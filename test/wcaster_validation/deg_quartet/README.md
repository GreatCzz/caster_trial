# Test 2 — 4 物种（必含 ref）：alignment_wcaster vs alignment_wtrial

**目的**：4 物种只有唯一四重奏且必含参考基因组 → wcaster 全局计分与 wtrial 的含-ref
计分应数学等价（验证加权 XXYY/scorePos 公式正确）。

**目录**：`deg_quartet/`（Test 2 of `test/wcaster_validation/`）

## 数据

`gen.py`（seed=20260904）生成 4 taxa（ref + a,b,c），2,000 bp 真实多样（无保守尾），权重非 1：
a≈0.9533, b≈0.9353, c≈0.8267, ref=1。

## 运行

```bash
./run.sh        # bin/alignment_wcaster 与 bin/alignment_wtrial, 同 fasta2ref
python3 verify.py
```

产物在 `results/`（两工具 .log/.tre）。

## 结果

| 判据 | 结果 |
|------|------|
| 物种权重 dump | 完全一致 |
| 全部 Score 行 | 16 行逐一相等，max diff = 0 |
| 最终树 RF | 0（均 `(((b,c),a),ref);`） |
