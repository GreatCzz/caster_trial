# wcaster 开发记录

系统发育物种树推断工具系列。`wcaster` = **CASTER + 序列相似度加权**（全局四重奏计分），分为两个加权粒度版本：
- `alignment_wcaster`（per-alignment 加权）
- `chunk_wcaster`（per-chunk 加权）

两者均已实现并通过退化测试（见 §5 与 `test/wcaster_validation/README.md`）。

---

## 1. 命名与文件

| 项目 | alignment_wcaster | chunk_wcaster |
|------|-------------------|---------------|
| 源文件 | `src/alignment_wcaster.hpp` | `src/chunk_wcaster.hpp` |
| 头文件保护 | `ALIGNMENT_WCASTER_HPP` | `CHUNK_WCASTER_HPP` |
| 命名空间 | `alignment_wcaster` | `chunk_wcaster` |
| 二进制 | `bin/alignment_wcaster` | `bin/chunk_wcaster` |
| 编译宏 | `-D ALIGNMENT_WCASTER` | `-D CHUNK_WCASTER` |
| 默认 chunk | 10000（同 CASTER） | 1000 |

> 历史：两个文件最初由 `caster.hpp` 复制改名；**alignment_wcaster 按 §3 改造为加权版后，
> chunk_wcaster 以 alignment_wcaster 为基复制（非再基于 caster）**，再按 §3.9 把权重改为 per-chunk。
> 文件名 / namespace / 编译后二进制名一致，命名方式与 wtrial 系列对应。

### 与现有工具的关系

| 维度 | CASTER | TRIAL 系(wtrial) | wcaster |
|------|--------|------------------|---------|
| 四重奏计分 | 全局对称（3 色任意取对） | 仅含 ref（ref 换入 slot0） | **全局对称（同 CASTER）** |
| 权重锚点 | 无 | 该 alignment 的参考物种 R | **该 alignment 的参考物种 R** |
| 参考物种优先入队 | 无 | 有 | **无（同 CASTER）** |
| 权重数组 | 无 | colorWeight/colorPairWeight | colorWeight/colorPairWeight |
| 输入 | fasta/phylip/列表 | fasta2ref | **fasta2ref** |
| 权重粒度 | — | alignment: per-alignment / chunk: per-chunk | alignment: per-alignment / chunk: per-chunk |

---

## 2. 数学基础

CASTER 整数 `XXYY(x0,x1,x2, y0,y1,y2)` 每一项为"某色组内取一对同类碱基 × 另两色各取一个单碱基"：

```
x0(x0-1)y1y2 + x1(x1-1)y2y0 + x2(x2-1)y0y1   ← X 型对在 color 0/1/2
+ y0(y0-1)x1x2 + y1(y1-1)x2x0 + y2(y2-1)x0x1  ← Y 型对在 color 0/1/2
```

加权后 `x_c(x_c-1)`（有序对）→ `2 × P(c)`，其中 `P(c)` 为该色组内碱基类 X 的**无序对**加权和
（`colorPairWeight` 采用无序约定，故 ×2 精确还原整数公式）。

> ⚠️ 与 wtrial 的关键区别：wtrial 的 ref 语义使对项只出现在 color1/2
> （color0 是 ref 换入的 slot），其 XXYY 只有 `x11,x22,y11,y22` 四个对参。
> **wcaster 全局计分里对项可出现在任意一个色组**，故签名需要 6 个对参
> （`x00,x11,x22,y00,y11,y22`）。命名沿用 wtrial 下标法。数学上必须如此才能退化为 CASTER。

### 二阶项约定

- `colorWeight[pos][color][nuc]` = 一阶加权和 `Σ w_i·c_i`
- `colorPairWeight[pos][color][nuc]` = 二阶同碱基无序对 `Σ_{i<j} w_i·w_j`（仅同碱基）
- 混合碱基对（嘌呤 A+G、嘧啶 C+T）用 wtrial 的 `pw()` 约定：
  - 单碱基对（A、G、C、T）：直接读 `cp[c][b]`
  - 嘌呤对 `pw(c,A,G) = cp[c][A] + cp[c][G] + cw[c][A]×cw[c][G]`
  - 嘧啶对 `pw(c,C,T) = cp[c][C] + cp[c][T] + cw[c][C]×cw[c][T]`
- `quadXXYY`/`quadPos`（NNI，每色取 1 个体）**无对项**，直接用 `colorWeight` 单和。

---

## 3. alignment_wcaster 实现（已完成）

源码改动均以醒目注释标记 `WCASTER-MOD ①..⑪`（chunk 版在复制 alignment 版后保留同名标记，
read 权重部分改用 `CHUNK-WCASTER-DEBUG` 包裹）。改动清单：

### 3.1 Attributes
`template<typename cnt_taxon_type = unsigned char>`，`cnt_t = cnt4_t = double`（删 `cnt_type` 参数，仿 wtrial）。
`cnt_taxon_t` 参与变体选择（bool/uchar/ushort 3 变体）。

### 3.2 SharedConstData::Element
加 `vector<cnt_t> speciesWeights;`（每 row 一个权重）。删除 `iReferenceTaxonId`
（全局计分不需要；ref 仅作 read 内权重锚点）。

### 3.3 Color 成员
`colorCnts` → 两个 double 数组：`colorWeight`（一阶 Σ w_i·c）、`colorPairWeight`（二阶同碱基无序对）。
构造器初始化两者为 `(data.nGenomePos)`。

### 3.4 elementSetOrClearTaxonColor（加权一阶/二阶维护）
```
w  = element.speciesWeights[row]
dw = (cnt_t)cnt * w
set  : colorPairWeight += colorWeight(旧) × dw ;  colorWeight += dw   // pair 先用旧 weight
clear: colorWeight -= dw ;  colorPairWeight -= colorWeight(新) × dw
```

### 3.5 XXYY（6 对参，命名仿 wtrial）
```
XXYY(x0,x1,x2, x00,x11,x22, y0,y1,y2, y00,y11,y22)
= 2·( x00·y1·y2 + x11·y2·y0 + x22·y0·y1 )
+ 2·( y00·x1·x2 + y11·x2·x0 + y22·x0·x1 )
```
`xNN` = color N 的 X-set 无序对加权和；`yN` = color N 的 Y-set 单和。

### 3.6 scorePos（全局对称，wtrial 写法）
输入 `(colorWeight, colorPairWeight, eqFreqs)`；预计算各色嘌呤/嘧啶/单碱基对
（`r00/r11/r22`、`y00/y11/y22`、`a*/g*/c*/t*` 含 color0）；9 次 XXYY 按 X/Y 碱基类传对参；
外层 `R2Y2` 组合式与 `pi` 系数同 CASTER。

### 3.7 quadXXYY / quadPos
公式不变；每色 1 个体无对项 → 参数 double，直接吃 `colorWeight` 单和。

### 3.8 elementScore / elementQuadripartitionScores
- `elementScore` → `scorePos(colorWeight[gPos], colorPairWeight[gPos], eqFreqs)`
- `elementQuadripartitionScores` → `quadPos(colorWeight[gPos], eqFreqs)`

### 3.9 read() — fasta2ref + 权重计算
输入为 **fasta2ref**（需每 alignment 的 ref 锚点）。权重：
```
alignment: 全序列非 gap 位点  Hamming → sim → w;  speciesWeights 赋给该 alignment 全部 chunk
chunk:     逐 chunk 在 [s,t) 全范围（含保守位点）算 → elements[iElementBegin+iChunk].speciesWeights[iRow]
w = (sim < 0.25) ? 0.0 : (sim - 0.25) / 0.75 ;  ref 权重 = 1.0
```
保留 bool/uchar/ushort 三档 `maxSpeciesman` 溢出检查；去掉 priorityTaxa/refColor 逻辑。

### 3.10 权重 dump（调试用）
- **alignment 版**：`WCASTER-MOD ⑩` 每 alignment 输出每物种单值权重（恒开，用于与 alignment_wtrial 对比；
  alignment_wtrial 侧 `WTRIAL-DEBUG-MOD` 同格式）。
- **chunk 版**：权重为 per-chunk（每物种每 chunk 不同），改为 `--dump-chunk-weights` flag **默认关闭**
  控制：输出 `chunk <idx> <species> <weight>`（按 chunk/物种名排序）。chunk_wtrial 亦加同 flag（`CHUNK-WTRIAL-DEBUG`）。

### 3.11 Driver
`DataClasses` → 3 变体（`<bool>/<uchar>/<ushort>`）；chunk 默认 10000（alignment）/ 1000（chunk）；
`programNames` 全名 "Alignment-Weighted CASTER" / "Chunk-Weighted CASTER"。

### 3.12 接入
`driver.cpp` 增加 `-D ALIGNMENT_WCASTER` / `-D CHUNK_WCASTER` 分支（include + namespace alias）。

---

## 4. 相关代码改动（wtrial 侧，供对比）

- `chunk_wtrial.hpp`：默认 chunk `10000`→`1000`；新增 `--dump-chunk-weights` flag（默认关闭）与 per-chunk dump。
- `alignment_wtrial.hpp`：`WTRIAL-DEBUG-MOD` 每 alignment 物种单值权重 dump（恒开）。

---

## 5. 测试（全部通过）

详见 `test/wcaster_validation/README.md`（含各子目录 README 与复现命令）：

| # | 目录 | 工具 | 结果 |
|---|------|------|------|
| 1 | `deg_weight1/` | alignment_wcaster vs CASTER（权重→1） | PASS：Score 接近、RF=0 |
| 2 | `deg_quartet/` | alignment_wcaster vs alignment_wtrial（4 物种） | PASS：Score 逐一相等、RF=0 |
| 3 | `chunk_deg_sim/` | chunk_wcaster vs chunk_wtrial（8 物种 100k per-chunk） | PASS：800 项权重 diff=0、RF=0 |
| 4 | `cat_chunk_compare/` | chunk_wcaster vs chunk_wtrial（cat，不 dump） | PASS：互比 + vs 真实树 RF 均 0 |
| 5 | `cat_weight_compare/` | alignment_wcaster vs alignment_wtrial（cat 权重） | PASS：10 物种权重 diff=0、RF=0 |

cat 参考真实树：`example/cat/true_species_tree.nwk`
（`cat_chunk_compare/make_true_tree.py` 从 CASTER cat 输出去 bootstrap 生成，供 cat 测试复用）。

## 6. 待办

- [ ] makefile 正式加入 `alignment_wcaster` / `chunk_wcaster` target（当前仅 driver.cpp 分支 + 手动 g++）
- [ ] CASTER.md / TRIAL.md 提及 wcaster 系列
- [ ] 大规模数据（如 570MB cat 的 alignment_wcaster vs caster 树级 sanity，可选）
