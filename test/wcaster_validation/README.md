# wcaster 测试验证

验证 `alignment_wcaster` 与 `chunk_wcaster` 的加权计分与权重计算正确性
（代码设计/修改见项目根 `wcaster.md`）。

## 测试索引

| # | 目录 | 工具 | 目的 | 核心判据 |
|---|------|------|------|----------|
| 1 | `deg_weight1/` | alignment_wcaster vs CASTER | 权重→1 退化 | Score 逐行接近、RF=0 |
| 2 | `deg_quartet/` | alignment_wcaster vs alignment_wtrial | 4 物种加权公式等价 | Score 逐一相等、RF=0 |
| 3 | `chunk_deg_sim/` | chunk_wcaster vs chunk_wtrial | per-chunk 权重一致（8 物种 100k） | (chunk,species) 权重 diff=0、RF=0 |
| 4 | `cat_chunk_compare/` | chunk_wcaster vs chunk_wtrial | cat 真实数据拓扑（不 dump） | RF=0（两工具互比 + vs 真实树） |
| 5 | `cat_weight_compare/` | alignment_wcaster vs alignment_wtrial | cat 真实数据权重一致 | 10 物种权重 diff=0、RF=0 |

每个子目录有自包含 `README.md`（数据/运行/结果）与 `results/` 产物。

## 汇总结果

| Test | 结果 |
|------|------|
| 1 | PASS：权重≈0.9992；16 行 Score 最大相对差 0.38%；RF=0 |
| 2 | PASS：权重一致；16 行 Score 逐一相等（diff=0）；RF=0 |
| 3 | PASS：800 项 (chunk,species) 权重集合一致、max\|diff\|=0；RF=0 |
| 4 | PASS：RF(chunk_wtrial,chunk_wcaster)=0；两者 vs true 均 =0 |
| 5 | PASS：10 物种权重 max\|diff\|=0；RF=0（bootstrap 差异源于计分范围不同，正常） |

参考真实树（Test 4/5 cat 复用）：`example/cat/true_species_tree.nwk`
（由 `cat_chunk_compare/make_true_tree.py` 从 CASTER cat 输出去 bootstrap 生成）。

## 复现命令

```bash
# Test 1（alignment 退化 vs CASTER）
cd test/wcaster_validation/deg_weight1        && python3 gen.py && ./run.sh && python3 verify.py
# Test 2（alignment 4 物种 vs wtrial）
cd test/wcaster_validation/deg_quartet        && python3 gen.py && ./run.sh && python3 verify.py
# Test 3（chunk 合成 8 物种 per-chunk 权重）
cd test/wcaster_validation/chunk_deg_sim      && python3 gen.py && ./run.sh
# Test 4（chunk cat 拓扑）
cd test/wcaster_validation/cat_chunk_compare  && python3 make_true_tree.py && ./run.sh
# Test 5（alignment cat 权重）
cd test/wcaster_validation/cat_weight_compare && ./run.sh
```

## 依赖

- `project/lib/tree_utils.py`（dendropy RF 距离）
- 先编译：`bin/caster`、`bin/alignment_wcaster`、`bin/alignment_wtrial`、`bin/chunk_wcaster`、`bin/chunk_wtrial`
