# caster_trial 项目

系统发育物种树推断工具系列。包含 ASTRAL、CASTER、SISTER、TRIAL、alignment_wtrial、chunk_wtrial、
以及 CASTER 的加权变体 alignment_wcaster、chunk_wcaster（见 `wcaster.md`）。

## 编译命令

```bash
make               # 编译全部工具（astral caster sister trial alignment_wtrial chunk_wtrial）
make trial         # 仅编译 TRIAL
make alignment_wtrial  # 仅编译 alignment_wtrial（per-alignment 权重）
make chunk_wtrial  # 仅编译 chunk_wtrial（per-chunk 权重）
```

- 编译器：g++ >= 13
- 标准：C++20，`-std=c++20 -march=native -Ofast`
- 所有工具共用 `src/driver.cpp`，通过预处理器宏选择编译：
  `-D ASTRAL`、`-D CASTER`、`-D SISTER`、`-D TRIAL`、`-D ALIGNMENT_WTRIAL`、`-D CHUNK_WTRIAL`
- wcaster 系列（`-D ALIGNMENT_WCASTER` / `-D CHUNK_WCASTER`）暂未加入 makefile，见 `wcaster.md`

## 项目结构

```
src/
├── driver.cpp                    # 主入口，通过 #ifdef 分派工具
├── driver.hpp                    # DRIVER concept 接口约束
├── common.hpp                    # 共享工具：LogInfo, InputParser, AnnotatedBinaryTree, Random
├── stepwise_colorable.hpp        # 接口概念：STEPWISE_COLORABLE, TAXON_ORDER_PRIORITIZING 等
├── optimization_algorithm.hpp    # 启发式搜索：TaxonOrderGenerator, subsample, placement, NNI, DP
├── placement_algorithm.hpp       # 逐步着色 placement 算法
├── nni_algorithm.hpp             # NNI 优化
├── constrained_dp_algorithm.hpp  # 约束动态规划树组装
├── quadripartition_support.hpp   # 四分区支撑度标注
├── alignment_utilities.hpp       # FASTA/Phylip 比对解析器
├── threadpool.hpp                # 多线程池
├── astral.hpp                    # ASTRAL 工具
├── caster.hpp                    # CASTER 工具
├── trial.hpp                     # TRIAL 工具（参考物种三角剖分）
├── alignment_wtrial.hpp          # alignment_wtrial 工具（per-alignment 加权）
├── chunk_wtrial.hpp              # chunk_wtrial 工具（per-chunk 加权）
├── alignment_wcaster.hpp         # alignment_wcaster 工具（CASTER + per-alignment 加权）
├── chunk_wcaster.hpp             # chunk_wcaster 工具（CASTER + per-chunk 加权）
├── sister.hpp                    # SISTER 工具
└── documentation.hpp             # DocumentationBase 基类
```

---

## TRIAL 算法

### 数据读取

1. 解析 `fasta2ref.txt`，每行 `<fasta路径> <参考物种名>`
2. 对每个 fasta 比对文件：两遍扫描（AP 统计频率做 chunking，AP2 填充核苷酸计数）
3. 产生 `Element` 数组，每个 element 包含一个 chunk 的有效位点和各物种的核苷酸计数

### 数据结构

```
SharedConstData
├── elements: vector<Element>
│   ├── Element.iGenomePosBegin    // 全局位点起始偏移
│   ├── Element.nPos               // 有效位点数
│   ├── Element.cnts[row][pos][nuc] // 物种×位点×4核苷酸的原始计数
│   ├── Element.taxon2row[taxId]   // 物种ID → row索引
│   ├── Element.eqFreqs[A/C/G/T]   // 平衡频率
│   └── Element.iReferenceTaxonId  // 该 element 的参考物种ID
├── nGenomePos                     // 全局位点总数
└── priorityTaxa                   // distinct 参考物种ID列表

Color (运行时)
├── colorWeight[pos][0..3][A/C/G/T]  // 加权和（非ref物种，weighted 版本）
├── colorPairWeight[pos][0..3][A/C/G/T] // 加权对和（非ref物种，weighted 版本）
├── refCnt[pos][A/C/G/T]            // 参考物种在各颜色的加权计数
└── refColor[pos]                   // 参考物种当前所在颜色组（-1=未设置）
```

### 评分：Tripartition（scorePos）

TRIAL 使用 3+1 颜色组方案：3 个非 ref 颜色组 + ref 所在颜色组。ref 的数据独立于 `colorWeight`，评分时通过 swap 将 ref 所在组换到位置 0：

```
elementScore(element):
  for pos in element:
    C = refColor[pos]           // ref 所在颜色
    cw = colorWeight[pos]       // 本地复制
    swap(cw[0], cw[C])          // ref 所在组 → 位置 0
    XXYY(refCnt, cw[0], cw[1], cw[2], pi)
    // ref 在位置 0，剩馀 3 组在 1,2,3（仅用 1,2 作为三色评分）
```

### 评分：Quadripartition（quadPos）

NNI 使用 4 颜色组评分，每色只取 1 个代表：

```
elementQuadripartitionScores(element):
  for pos:
    C = refColor[pos]
    cw = colorWeight[pos]
    cw[C] = refCnt[pos]         // ref 替代该颜色组
    quadPos(cw, pi)             // 四组均等，原始 CASTER 公式
```

### 参考物种优先入队

通过 `TAXON_ORDER_PRIORITIZING` 概念和 `TaxonOrderGenerator` 自动将 `priorityTaxa` 中的所有物种排到随机序列最前面，保证 tree building 过程中 ref 最先加入树。

### Element::cnts 类型回退链

```
DataClasses = variant<Color<StepwiseColorDefaultAttributes<bool>>::SD,
                      Color<StepwiseColorDefaultAttributes<uchar>>::SD,
                      Color<StepwiseColorDefaultAttributes<ushort>>::SD>

getStepwiseColorSharedConstData():
  try read<变体0>()  → throw（bool不够）? → catch
  try read<变体1>()  → throw（uchar不够）? → catch
  return read<变体2>()
```

3 个变体而非原始 CASTER 的 5 个（`cnt_t` 固定为 `double`，无需组合 `cnt_type` 变体）。

---

## alignment_wtrial 算法（per-alignment 加权）

### 权重计算

对每个 alignment，计算所有非 ref 物种的汉明距离→相似度→线性权重：

```
similarity = 1 - Σ[ref≠species且非gap] / Σ[非gap]
weight = 0                          (similarity < 0.25)
       = (similarity - 0.25) / 0.75 (similarity ≥ 0.25)
ref权重 = 1.0

所有 chunk 共享同一权重 → speciesWeights[row] = weight
```

### 评分：colorPairWeight

TRIAL 的 `y2*(y2-1)` 组合数公式对加权小数计数不成立。weighted 版本引入独立的加权对和数组，在评分中替代组合项：

```
colorWeight[pos][color][nuc]       // Σ(w_i × c_i)
colorPairWeight[pos][color][nuc]   // Σ_{i<j} (w_i×c_i) × (w_j×c_j)
```

**增量更新**（`elementSetOrClearTaxonColor`）：

| 操作 | 公式 |
|------|------|
| 加入物种（权重 w，计数 c） | `pairWgt += oldWeight × w×c`；`weight += w×c` |
| 移除物种 | `weight -= w×c`；`pairWgt -= newWeight × w×c` |

注意 set 时先更新 `pairWeight` 再 `weight`，clear 时先更新 `weight` 再 `pairWeight`。

**XXYY 签名**（+4 个 pair-weight 参数）：

```
XXYY(xR, x0, x1, x2, yR, y0, y1, y2, x11, x22, y11, y22)
```

- `x11 = pairWgt[1][A] + pairWgt[1][G] + wgt[1][A]×wgt[1][G]`（颜色 1 嘌呤对）
- `x22 = pairWgt[2][A] + pairWgt[2][G] + wgt[2][A]×wgt[2][G]`
- `y11 = pairWgt[1][C] + pairWgt[1][T] + wgt[1][C]×wgt[1][T]`
- `y22 = pairWgt[2][C] + pairWgt[2][T] + wgt[2][C]×wgt[2][T]`

同核苷酸对（如 AA）直接从 `pairWgt[color][nuc]` 读取，无需 `pw()` 函数。

**公式体**：

```
return xR×x0×y1×y2×2 + yR×y0×x1×x2×2       // 跨组单取
     + xR×(xR-1)×y1×y2 + yR×(yR-1)×x1×x2    // ref 自对
     + xR×x1×y22×2 + yR×y1×x22×2             // ref + col1, col2自对
     + xR×x2×y11×2 + yR×y2×x11×2             // ref + col2, col1自对
```

### 数据结构（weighted 特有）

```
Element.speciesWeights[row]   // 物种×1，该 element 的权重（per-alignment 版本全 chunk 相等）
Color.colorWeight[pos][0..3][A/C/G/T]
Color.colorPairWeight[pos][0..3][A/C/G/T]
```

### 验证

| 阶段 | 数据集 | 判据 | 结果 |
|------|------|------|------|
| Phase 1 | 4 物种 2000bp，权重 ≈ 1.0 | RF=0，bootstrap 逐位相等 | PASS：quadPos 公式等价 |
| Phase 2 | Python XXYY 独立实现 vs C++ 同一输入 | 数学一致 | PASS：6 组已知输入输出验证 |
| Phase 3 | score_validation 数据集 | `alignment_wtrial_score / trial_score = weight_product` | PASS：值等于 0.601 |

---

## chunk_wtrial 算法（per-chunk 加权）

与 alignment_wtrial 唯一区别在**权重计算粒度的变化**。其余完全相同（评分、数据结构、类型系统）。

### 权重计算

```
while (AP2.nextSeq()):
  for chunk in 0..nChunk-1:
    // 在 chunk 全范围 [s, t) 的所有位点上计算汉明距离（含保守位点）
    s = chunk * nSites / nChunk, t = (chunk+1) * nSites / nChunk
    for iPos in [s, t):
      if (refSeq[iPos]=='-' || seq[iPos]=='-') continue
      nonGap++; if (refSeq != seq) hamming++
    nonGap == 0 → weight = 0    // 全 gap 处理
    否则 weight = (sim - 0.25) / 0.75

    // 权重直接赋值给每个 element
    sharedConstData.elements[iElementBegin + chunk].speciesWeights[row] = weight
```

**注意**：相似度统计遍历 **chunk 全范围所有位点**（含保守位点），而非仅 `sites[chunk]` 的多样性位点。`cnts` 填充仍只使用多样性位点（`sites[chunk]`），这是 CASTER 的评分信号设计。此 bug（仅用多样性位点统计相似度）于 2026-08-22 修复。

### 与 alignment_wtrial 的对比

| 特性 | alignment_wtrial | chunk_wtrial |
|------|--------|-------------|
| 权重粒度 | per-alignment | **per-chunk** |
| 源文件 | `alignment_wtrial.hpp` | `chunk_wtrial.hpp` |
| 编译宏 | `-D ALIGNMENT_WTRIAL` | `-D CHUNK_WTRIAL` |
| 二进制 | `bin/alignment_wtrial` | `bin/chunk_wtrial` |
| 命名空间 | `alignment_wtrial` | `chunk_wtrial` |
| 默认 chunk | 10000 | 1000（`--chunk` 可配） |
| 评分 | colorPairWeight 精确 | colorPairWeight 精确 |
| 相似度统计范围 | 全序列（无 bug） | 全 chunk 范围（修复前仅多样性位点） |
| 数学正确性 | ✓ | ✓ |

### chunk 大小对精度的影响（10 物种猫科数据，2026-08-22 修复后）

| chunk | RF vs CASTER | 结论 |
|-------|-------------|------|
| 1000 | 0 | ✓ 与 CASTER 一致 |
| 2000 | 0 | ✓ 与 CASTER 一致 |
| 5000 | 0 | ✓ 与 CASTER 一致 |
| 10000 | 0 | ✓ 与 CASTER 一致 |
| 50000 | 0 | ✓ 与 CASTER 一致 |

修复前 chunk≤2000 的 RF=8~10、chunk=50000 的 RF=2，修复后全部归零。
bootstrap 随 chunk 增大略增（Leopardus/Caracal 分支 92.6→96.3），说明更大 chunk 权重更稳。

默认 chunk 现为 1000（与 chunk_wcaster 一致）；需要更粗/更稳权重时可调大 `--chunk`（如 5000/10000）。

---

## 类型系统

### cnt 三层体系

| 类型 | 作用 | weighted 值 |
|------|------|---------|
| `cnt_taxon_t` | `Element::cnts` 原始计数 | `bool`/`uchar`/`ushort`（3 变体） |
| `cnt_t` | `colorWeight`/`colorPairWeight` | `double`（固定） |
| `cnt4_t` | XXYY 乘法中间结果 | `double`（固定） |

### 模板参数

```cpp
template<typename cnt_taxon_type>  // bool / uchar / ushort
struct StepwiseColorDefaultAttributes {
    using cnt_taxon_t = cnt_taxon_type;
    using cnt_t  = double;
    using cnt4_t = double;
    using index_t = long long;
    using score_t = double;
};
```

只有 `cnt_taxon_t` 参与变体选择（3 个 DataClass），`cnt_t` 和 `cnt4_t` 固定为 `double`。

---

## 测试数据

```
test/
├── pairweight_validation/      # 评分验证
│   ├── phase1_weight1/         # RF=0 + bootstrap 逐位验证（quadPos 正确）
│   └── phase2_manual/          # Python XXYY 独立实现 + 单元测试
├── score_validation/           # 数字验证（weighted/TRIAL = 权重乘积）
├── debug_mapping/              # 集群 -a mapping 诊断
└── wcaster_validation/         # wcaster 退化测试（deg_weight1/deg_quartet/cat_weight_compare/
                               #   chunk_deg_sim/cat_chunk_compare），见 wcaster.md
example/
├── cat/                        # 猫科全基因组（570 MB, 10 Felidae）
│   ├── test_full.fasta
│   ├── fasta2ref.txt
│   ├── true_species_tree.nwk   # 参考真实树（去 bootstrap，cat 测试复用）
│   └── results_*/              # CASTER / TRIAL / weighted 对比
└── test_small/
```

## 工具库

- `../lib/tree_utils.py` — 基于 dendropy 的 RF 距离、树清洗、剪枝
  使用 `rf_distance(t1, t2)` 验证拓扑一致性（RF=0 为相同树）
