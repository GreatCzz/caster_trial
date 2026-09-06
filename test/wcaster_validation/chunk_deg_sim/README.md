# Test 3 — 8 物种 100k 合成数据：chunk_wcaster vs chunk_wtrial per-chunk 权重

**目的**：chunk 版两工具（同一权重公式/输入）per-chunk 权重逐项一致，且拓扑一致。
**目录**：`chunk_deg_sim/`（Test 3 of `test/wcaster_validation/`）

## 数据

`gen.py`（seed=20260906）：ref + a..g 共 8 物种，100,000 bp，`--chunk 1000`（100 chunk）。
各物种按 chunk 变化的替换率模拟 → per-chunk 权重逐 chunk 变化、各物种不同
（均值 0.62~0.92，范围约 0.42~0.97），非全 1 值。

## 运行

```bash
./run.sh        # chunk_wtrial 与 chunk_wcaster, --chunk 1000, 均 --dump-chunk-weights
```

产物在 `results/`。`compare_weights.py` 解析两 log 的 "Chunk weights for alignment file:" 块，
按 (chunk,species) 比对；`compare_topology.py` 算两树 RF。

## 结果

| 判据 | 结果 |
|------|------|
| (chunk,species) 权重 | 800 项集合完全一致，max\|diff\| = 0 |
| 拓扑 RF | 0（两工具一致） |
