# Test 4 — cat 真实数据：chunk_wcaster vs chunk_wtrial 拓扑（不输出权重）

**目的**：cat 真实数据下两 chunk 工具拓扑一致，且均与参考真实树一致。
**目录**：`cat_chunk_compare/`（Test 4 of `test/wcaster_validation/`）

## 参考真实树

`example/cat/true_species_tree.nwk`：由 `make_true_tree.py` 复用
`results_20260722_1/caster.tre` 去 bootstrap/枝长生成（真实物种名），后续 cat 测试复用，无需重跑 CASTER。

```bash
python3 make_true_tree.py   # 生成/刷新 example/cat/true_species_tree.nwk
```

## 运行

```bash
./run.sh        # chunk_wtrial 与 chunk_wcaster, --chunk 1000（不传 dump flag）
```

产物在 `results/`。`compare_topology.py` 计算：
- RF(chunk_wtrial, chunk_wcaster)
- RF(chunk_wtrial, true)、RF(chunk_wcaster, true)

## 数据

`example/cat/test_full.fasta`（570 MB, 10 Felidae），ref=Felis_catus，root=Canis_lupus_familiaris。

## 结果

| RF | 值 |
|----|----|
| chunk_wtrial vs chunk_wcaster | 0 |
| chunk_wtrial vs true | 0 |
| chunk_wcaster vs true | 0 |
