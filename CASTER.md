# caster_trial 项目

系统发育物种树推断工具系列。CASTER 是核心，其余工具在其基础上扩展。

## 编译命令

```bash
make caster       # 仅编译 CASTER
make              # 编译全部工具
```

- 编译器：g++ >= 13
- 标准：C++20，`-std=c++20 -march=native -Ofast`

## 项目结构

```
src/
├── driver.cpp                    # 主入口，#ifdef 分派工具
├── common.hpp                    # LogInfo, InputParser, AnnotatedBinaryTree, Random, ChangeLog
├── stepwise_colorable.hpp        # Concept: STEPWISE_COLORABLE, QUADRIPARTITION_COLORABLE
├── optimization_algorithm.hpp    # 启发式搜索
├── placement_algorithm.hpp       # 逐步着色放置算法
├── nni_algorithm.hpp             # NNI 优化
├── constrained_dp_algorithm.hpp  # 约束动态规划
├── alignment_utilities.hpp       # FASTA/Phylip 解析
├── threadpool.hpp                # 多线程池
├── caster.hpp                    # ★ CASTER-site
├── trial.hpp                     # TRIAL（参考三角剖分）
├── alignment_wtrial.hpp          # alignment_wtrial（per-alignment 加权）
├── chunk_wtrial.hpp              # chunk_wtrial（per-chunk 加权）
└── ...
```

TRIAL、alignment_wtrial、chunk_wtrial 的详细说明见 `TRIAL.md`。

---

## CASTER 算法

### 输入

一个多序列比对文件（FASTA 或 Phylip），每物种一个或多个单倍体。

### 数据读取（`DriverHelper::read()`）

两遍扫描（`AlignmentParser AP` + `AP2`）：

**第一遍**（`AP.nextSeq()`）：统计每个位点的 A/C/G/T 频率 → 确定有效位点（purine ≥ 2 且 pyrimidine ≥ 2 的位点）→ 分割为 chunk → 计算每个 chunk 的平衡频率（eqFreqs）

**第二遍**（`AP2.nextSeq()`）：只读有效位点，填入 `Element::cnts[row][pos][nuc]`（每物种·每位点·每核苷酸的原始计数，0 或 1）

**chunking**：按 `--chunk` 参数（默认 10000）将序列按位点数切为多个 Element。每个 Element 是评分计算的基本单位。

**类型选择**：5 个变体（`bool/uchar`、`uchar/uchar`、`bool/ushort`、`uchar/ushort`、`ushort/ushort`），自动选择最小的够用的整数类型存 `cnts` 和 `colorCnts`。

### 数据结构

```
SharedConstData::Element
├── iGenomePosBegin      // 全局位点起始偏移
├── nPos                 // 有效位点数
├── cnts[row][pos][nuc]  // 原始计数（cnt_taxon_t：0或1）
├── taxon2row[taxId]     // 物种ID → row索引
└── eqFreqs[A/C/G/T]     // 该chunk的平衡频率

Color (运行时)
├── sharedConstData&      // 不可变数据引用
├── colorCnts[pos][0..3][A/C/G/T]  // 4 颜色组的核苷酸计数（cnt_t）
└── refColor[pos]          // TRIAL/weighted 版本特有
```

### 颜色分配（placement 算法）

逐步放置（Recursive Placement）或两步放置（TwoStep Placement）算法依次将物种加入树。每一步对一个内部节点进行后序遍历，将物种分配为 3 个颜色组：

```
colorSubtreeSet(node):
  if node is leaf:
    setColor(leaf, 2)     // 所有叶子初始为 color0，移动到 color2
  
  else:
    rightSubtreeProcess(node.right)
    clearSetColor(right, 2, 0)   // 右子树：2→0
    leftSubtreeProcess(node.left)
    clearSetColor(right, 0, 1)   // 右子树：0→1
    
    computeScore()  // ★ 此时 3 色分配为：
                    //  color0 = 当前节点外侧
                    //  color1 = 右子树
                    //  color2 = 左子树
    
    clearSetColor(right, 1, 2)   // 右子树：1→2
    // 合并左右子树

  clearSetColor(newTaxon, 0, 1)
  computeScore()  // NEW_NODE
```

**3 个颜色组的语义**：

| 颜色 | 含义 |
|------|------|
| 0 | 当前节点外侧（与左右子树并列的分支） |
| 1 | 右子树 |
| 2 | 左子树 |

颜色 3 仅在四分区（quadripartition）评分中使用（NNI），placement 不设置。

### 评分：Tripartition（scorePos）

对每个位点，3 个颜色组的核苷酸统计形成三色四叶树模式。公式为原始 CASTER 的 `XXYY`：

```
elementScore(element):
  for pos in element:
    c = colorCnts[pos]            // c[0..2][A/C/G/T] — 3 个颜色组的核苷酸计数
    res += scorePos(c, eqFreqs)

scorePos(c, pi):
  a0 = c[0][A], c0 = c[0][C], g0 = c[0][G], t0 = c[0][T]   // 颜色0
  a1 = c[1][A], ...                                           // 颜色1
  a2 = c[2][A], ...                                           // 颜色2

  r0 = a0 + g0, y0 = c0 + t0    // 嘌呤/嘧啶和
  r1 = a1 + g1, y1 = c1 + t1
  r2 = a2 + g2, y2 = c2 + t2

  // 9 个 XXYY 调用（所有核苷酸组合）
  rryy = XXYY(r0, r1, r2, y0, y1, y2)
  aayy = XXYY(a0, a1, a2, y0, y1, y2)
  ggyy = XXYY(g0, g1, g2, y0, y1, y2)
  rrcc = XXYY(r0, r1, r2, c0, c1, c2)
  rrtt = XXYY(r0, r1, r2, t0, t1, t2)
  aacc = XXYY(a0, a1, a2, c0, c1, c2)
  aatt = XXYY(a0, a1, a2, t0, t1, t2)
  ggcc = XXYY(g0, g1, g2, c0, c1, c2)
  ggtt = XXYY(g0, g1, g2, t0, t1, t2)

  return rryy×R2×Y2 -(aayy+ggyy)×R²×Y² -(rrcc+rrtt)×R2×Y² +(aacc+aatt+ggcc+ggtt)×R²×Y²
```

### XXYY 公式

```cpp
XXYY(x0, x1, x2, y0, y1, y2) =
    x0*(x0-1)*y1*y2 + x1*(x1-1)*y2*y0 + x2*(x2-1)*y0*y1 +
    y0*(y0-1)*x1*x2 + y1*(y1-1)*x2*x0 + y2*(y2-1)*x0*x1
```

**语义**：从 3 个颜色组中各取若干嘌呤/嘧啶类型个体，计算所有满足"第一对属于 X 型核苷酸、第二对属于 Y 型"的四叶树组合数之总和。

公式的对称性：
- 对任意两个位置的交换是非对称的（`x0 vs x1` 不相等）
- 对循环排列（0→1→2→0）是对称的——所有项均包含全 3 个颜色组

### 评分：Quadripartition（quadPos）

NNI 优化的评分。每条内部边的 4 个分支（色 0/1/2/3）各自作为独立的颜色组。使用 `quadXXYY`：

```cpp
quadXXYY(x0,x1,x2,x3, y0,y1,y2,y3) = x0*x1*y2*y3 + y0*y1*x2*x3
```

每色只取 1 个个体，无 pair 项。四分区得分用于：
1. NNI 算法判断 "围绕颜色 0 节点的三种替代拓扑中哪个得分最高"
2. Bootstrap 标注

### 启发式搜索（optimization_algorithm）

```
heuristSearch(data, nTaxa, r, s, nThreads):
  Random random
  Color color(data)
  树 = subsampleProcedure2(color, r rounds)
  do:
    树 = subsampleProcedure2(color, s rounds)
  while (newScore > oldScore + EPSILON)
  return 树
```

**subsampleProcedure2**：随机子采样。从全部 element 中随机抽取 min(`subsample-min`, total) 个 element，递归 placement 构建树 → NNI 优化 → 约束 DP 聚合 → 输出最佳树。

**r** 和 **s** 参数控制迭代轮数：首轮 r 次，后续 s 次。

---

### CASTER-site 得分含义

`scorePos` 计算的是 CASTER-site 的**位点模式得分**（site pattern score），而非逐棵四叶树的计数。该得分对所有四叶树是**线性求和**的——3 颜色组的分配决定了对树中每个内部节点的支撑评价，所有内部节点的得分总和即为整棵树的 CASTER-site 得分。

---

## 类型体系

| 类型 | 角色 | CASTER 中的值 |
|------|------|-------------|
| `cnt_taxon_t` | `Element::cnts` 原始计数 | `bool`/`uchar`/`ushort`（每物种 0~多个样本） |
| `cnt_t` | `colorCnts` 颜色组累加 | `uchar`/`ushort`（多个物种累加） |
| `cnt4_t` | XXYY 乘法 | `uint`/`ull`（乘法溢出保护） |
| `score_t` | 最终得分 | `double` |
| `index_t` | 位点索引 | `long long` |

5 变体 DataClasses 通过 `try-catch` 自动回退选择最优内存类型。
