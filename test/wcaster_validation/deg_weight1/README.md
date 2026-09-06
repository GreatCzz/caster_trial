# Test 1 — 权重→1 退化：alignment_wcaster vs CASTER

**目的**：权重≈1 时，wcaster 加权计分应退化为 CASTER 整数计分（验证加权 XXYY/scorePos）。

**目录**：`deg_weight1/`（Test 1 of `test/wcaster_validation/`）

## 数据

`gen.py`（seed=20260904）生成 7 taxa（ref + s1..s6），每序列 1,000,800 bp：
- 前 800 bp：多样性块（随机 ACGT，信息位点）
- 后 1,000,000 bp：全 A 保守尾（长尾把相似度推向 ≈1，且不与计分位点同 chunk）

权重：s1..s6 ≈ 0.9992，ref = 1。

## 运行

```bash
./run.sh        # bin/caster 与 bin/alignment_wcaster, --chunk 800 -t 1 --initial-round 4 --subsequent-round 2
python3 verify.py
```

`run.sh` 产物在 `results/`（两工具 .log/.tre）。

## 结果

| 判据 | 结果 |
|------|------|
| 物种权重 ≈ 1 | PASS（s1..s6 = 0.99918~0.99923） |
| 逐行 Score 相对差 | PASS（16 行最大相对差 0.38%） |
| 最终树 RF vs CASTER | 0（`(((((s4,s1),s6),s3),(s2,s5)),ref);`） |

微小相对差来自权重≈0.9992 而非严格 1，符合预期。
