# debug_mapping — alignment_wtrial `-a mapping` crash

## 问题

alignment_wtrial 在集群上传入 `-a mapping.txt` 参数时崩溃（Illegal instruction），而 CASTER 和 TRIAL 均正常。

## 本地验证结果（2026-07-23）

| Bug | 本地状态 | 说明 |
|------|---------|------|
| alignment_wtrial + `-a` 崩溃 | **不重现** | 本地 `(((6,5),9),4),(8,10)),7),3),2),1)` 正常输出，标签映射正确 |
| alignment_wtrial 不带 `-a`，标签为 `genome_X` | 预期行为 | CASTER/TRIAL 同样需要 `-a` 做标签映射 |
| alignment_wtrial 不带 `-a`，末尾多出 `,1)` 叶子 | 预期行为 | fa2ref 中 ref="1" 不在 fasta（fasta 中是 `genome_1`），需 `-a` 映射 |

**集群崩溃原因推测**：集群 CPU 指令集与本地不同。可能的原因：
1. 集群编译时使用的 `-march=native` 产生了本地不支持的指令（如 AVX-512）
2. 集群上的 `alignment_wtrial` 二进制来自本地，通过 scp 传输后指令集不兼容

**解决方案**：在集群上重新编译，去掉 `-march=native`：
```bash
g++ -std=c++20 -O3 -D ALIGNMENT_WTRIAL src/driver.cpp -o bin/alignment_wtrial
```
或使用兼容性更好的 `-march=x86-64-v2`。

## 文件说明

| 文件 | 用途 |
|------|------|
| `rep_1_ref_genome_1_100pct.fasta` | r5 rep_1 对齐文件（ref=genome_1, 10 species） |
| `mapping.txt` | 标签映射 `genome_1 → 1` ~ `genome_10 → 10` |
| `fa2ref_mr_1.txt` | 单 ref (k=1) 输入 |
| `fa2ref_mr_3.txt` | 多 ref (k=3) 输入 |
| `truth_s_tree.nwk` | 真物种树（含外群 0） |

## 对比

### CASTER（正常）

```bash
caster -i aln.fasta -a mapping.txt -t 8 --root 10 -o out.tre
# 产出标签: 1,2,3,...,10
```

### TRIAL（正常）

```bash
trial -i fa2ref.txt -a mapping.txt -t 8 -o out.tre
# 产出标签: 1,2,3,...,10
```

### alignment_wtrial（集群崩溃，本地正常）

```bash
alignment_wtrial -i fa2ref.txt -a mapping.txt -t 8 -o out.tre
# 集群: Illegal instruction (core dumped)
# 本地: 正常输出，标签映射正确
```

### alignment_wtrial 不带 -a 的行为

```bash
alignment_wtrial -i fa2ref.txt -t 8 -o out.tre
# 产出的树：
# 1. 标签为 genome_1,genome_2,...（未映射，需 -a）
# 2. fa2ref 第二列的 ref ID 被当作额外叶子节点 → 末尾多出 ,1),2),3)
# 3. 所有枝长为 0
# 以上均为预期行为：fa2ref 中的 ref 名 "1" 不匹配 fasta 中的 "genome_1"
```

## 期待行为

alignment_wtrial 应接受 `-a mapping.txt` 参数，输出含数字标签（1,2,...）且不含多余 ref ID 叶子的树。（本地已验证正常）
