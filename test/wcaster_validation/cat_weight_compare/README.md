# Test 5 — cat 数据物种权重对比：alignment_wtrial vs alignment_wcaster

**目录**：`cat_weight_compare/`（Test 5 of `test/wcaster_validation/`；总索引见上级 `README.md`）

**目的**：在真实 cat 数据集上对比 `alignment_wtrial` 与 `alignment_wcaster` 的 per-alignment 物种权重计算是否一致。
两者权重公式相同（以各 alignment 的参考物种为锚点，Hamming 距离 → similarity → 线性权重），
wcaster 是 CASTER 派生（全局计分），wtrial 是 TRIAL 派生（ref 限制计分），但权重计算逻辑应完全一致。

**数据**：`example/cat/test_full.fasta`（570 MB，10 物种，Felidae），ref = Felis_catus，root = Canis_lupus_familiaris。

## 文件

| 文件 | 作用 |
|------|------|
| `run.sh` | 运行两工具并调用对比（数据经 `../../../example/cat/fasta2ref.txt`） |
| `compare_weights.py` | 解析两 log 的 "Species weights" dump，逐物种 diff（判据：物种集合相同且 max\|Δ\|<1e-9） |
| `results/alignment_wtrial.{log,tre}` | alignment_wtrial 输出 |
| `results/alignment_wcaster.{log,tre}` | alignment_wcaster 输出 |
| `results/weight_compare.txt` | 逐物种权重对照表 |

> 权重 dump 说明：wtrial 的 `alignment_wtrial.hpp` 含 WTRIAL-DEBUG-MOD 调试输出，
> wcaster 的 `alignment_wcaster.hpp` 含 WCASTER-MOD ⑩ 调试输出，格式一致：
> `Species weights for alignment file: ...` + 每物种 `  <name> = <weight>`。

## 运行命令

```bash
cd test/wcaster_validation/cat_weight_compare && ./run.sh
```

参数：`-t 8 --initial-round 4 --subsequent-round 2 --root Canis_lupus_familiaris`（同 example/cat 历史 run）。

## 测试结果（2026-09-06）

权重 dump **10 物种逐一完全一致**（见 `results/weight_compare.txt`），max |diff| = 0：

| species | alignment_wtrial | alignment_wcaster | diff |
|---------|------------------|-------------------|------|
| Acinonyx_jubatus | 0.9737860000 | 0.9737860000 | 0 |
| Canis_lupus_familiaris | 0.7636140000 | 0.7636140000 | 0 |
| Caracal_caracal | 0.9725480000 | 0.9725480000 | 0 |
| Felis_catus | 1.0000000000 | 1.0000000000 | 0 |
| Leopardus_geoffroyi | 0.9705930000 | 0.9705930000 | 0 |
| Lynx_canadensis | 0.9754110000 | 0.9754110000 | 0 |
| Otocolobus_manul | 0.9759040000 | 0.9759040000 | 0 |
| Panthera_tigris | 0.9666000000 | 0.9666000000 | 0 |
| Prionailurus_bengalensis | 0.9763370000 | 0.9763370000 | 0 |
| Suricata_suricatta | 0.8275160000 | 0.8275160000 | 0 |

**VERDICT: PASS**（物种集合相等，max\|diff\| = 0）

### 树 sanity（RF 距离）

`rf_distance(alignment_wtrial.tre, alignment_wcaster.tre) = 0` —— 两工具最终拓扑一致。

两树仅在 bootstrap 上有正常差异（wtrial 95.8 vs wcaster 99.3 @ Leopardus/Caracal 分支），
源于计分范围不同（wtrial ref 限制 vs wcaster 全局），与权重无关。
