# chunk_wtrial — chunk 大小对 RF 精度的影响

> ## ⚠️ 本结果已失效（2026-08-01 运行，标记于 2026-08-22）
>
> 该批次测试的权重相似度统计存在 **bug**：相似度（汉明距离）只基于
> **多样性位点**（`sites[iChunk]`，即 `freq[A]+freq[G]≥2 && freq[C]+freq[T]≥2`
> 的位点）计算，排除了大量保守的、比对良好的区段，导致相似度被严重低估、
> 差异被严重高估。
>
> 该 bug 已于 **2026-08-22** 在 `chunk_wtrial.hpp` 中修复：
> 权重计算改为遍历 chunk 全范围 `[s, t)` 的所有位点（含保守位点）。
>
> 正确结果请参见 `../chunk_test_results_2026-08-22.md`（重新运行）。

## 测试目的

评估 `--chunk` 参数（控制每个 element 的有效位点数）对 chunk_wtrial 推断树拓扑精度的影响。chunk 越小 → element 越细 → per-chunk 权重噪声越大；chunk 越大 → element 越粗 → 权重估算越稳。

## 数据集

| 属性 | 值 |
|------|------|
| 文件 | `example/cat/test_full.fasta` |
| 大小 | 570 MB |
| 物种 | 10 个（Felidae 猫科，含 Canis_lupus_familiaris 和 Suricata_suricatta 作为外群） |
| 位点数 | ~500,000 bp |
| 参考物种 | Felis_catus |
| 置根 | Canis_lupus_familiaris（`--root Canis_lupus_familiaris`） |

## 运行参数

所有测试使用相同参数，仅 `--chunk` 不同：

```bash
chunk_wtrial -i fasta2ref.txt --chunk <N> -t 4 --initial-round 4 --subsequent-round 2 \
  --root Canis_lupus_familiaris --log chunk_<N>.log -o chunk_<N>.tre
```

基准树由 CASTER 生成（`-D CASTER`，默认 chunk=10000，无加权）：

```
((((((((Prionailurus_bengalensis,Otocolobus_manul)100,Felis_catus)100,
Lynx_canadensis)100,Acinonyx_jubatus)100,Leopardus_geoffroyi)98.9,
Caracal_caracal)100,Panthera_tigris)100,Suricata_suricatta),Canis_lupus_familiaris);
```

## 结果

### RF 距离汇总

| chunk | RF vs CASTER | 元素数 | 每元素平均有效位点 | 结论 |
|-------|:-----------:|--------|-------------------|------|
| 1000 | **10** | ~500 | ~1,000 | 每 chunk 仅 ~1000 位点，其中有效位点更少，权重噪声极大，树拓扑完全偏离 |
| 2000 | **8** | ~250 | ~2,000 | 仍不稳定，Felidae 内部分支混乱 |
| **5000** | **0** | ~100 | ~5,000 | ✓ 与 CASTER 完全一致 |
| **10000** | **0** | ~50 | ~10,000 | ✓ 与 CASTER 完全一致 |
| 50000 | **2** | ~10 | ~50,000 | 退化为 per-alignment（element 数太少），Suricata 和 Panthera 姐妹群错位 |

### 每 chunk 输出树

| chunk | RF | tree |
|-------|-----|------|
| 1000 | 10 | ((((((((Caracal_caracal,Panthera_tigris)83.6,Leopardus_geoffroyi)100,Acinonyx_jubatus)99,Lynx_canadensis)100,(Prionailurus_bengalensis,Otocolobus_manul)100)86.4,Felis_catus)100,Suricata_suricatta),Canis_lupus_familiaris); |
| 2000 | 8 | (((((((Panthera_tigris,Caracal_caracal)92.5,Leopardus_geoffroyi)100,Acinonyx_jubatus)99.1,Lynx_canadensis)94.9,((Prionailurus_bengalensis,Otocolobus_manul)100,Felis_catus)99.5)100,Suricata_suricatta),Canis_lupus_familiaris); |
| 5000 | 0 | (((((((((Otocolobus_manul,Prionailurus_bengalensis)100,Felis_catus)100,Lynx_canadensis)99.7,Acinonyx_jubatus)100,Leopardus_geoffroyi)97.2,Caracal_caracal)99.6,Panthera_tigris)100,Suricata_suricatta),Canis_lupus_familiaris); |
| 10000 | 0 | (((((((((Otocolobus_manul,Prionailurus_bengalensis)100,Felis_catus)100,Lynx_canadensis)100,Acinonyx_jubatus)100,Leopardus_geoffroyi)99.3,Caracal_caracal)99.7,Panthera_tigris)88.6,Suricata_suricatta),Canis_lupus_familiaris); |
| 50000 | 2 | ((((((((Prionailurus_bengalensis,Otocolobus_manul)100,Felis_catus)100,Lynx_canadensis)100,Acinonyx_jubatus)100,Leopardus_geoffroyi)99.6,Caracal_caracal)95.1,(Suricata_suricatta,Panthera_tigris)69.1),Canis_lupus_familiaris); |

## 性能对比

| chunk | 耗时 | 内存 | element 数 |
|-------|------|------|-----------|
| 1000 | ~28s | ~615 MB | ~500 |
| 2000 | ~25s | ~615 MB | ~250 |
| 5000 | ~23s | ~615 MB | ~100 |
| 10000 | ~22s | ~615 MB | ~50 |
| 50000 | ~20s | ~615 MB | ~10 |

chunk 数与耗时成正比（更多 element 需要更多次 scorePos 调用），但差异不大（28s vs 20s）。

## 日志内容

每份日志（`.log` 文件）包含：
- 输入参数回显
- 数据解析统计（#Taxa、#Elements、#Threads）
- 每轮 subsample 的树构建过程（Recursive placement → NNI → Score）
- 最终树和 bootstrap 标注
- 时间统计

## 结论

1. **推荐默认 chunk = 10000**：RF = 0，所有 benchmark 均与 CASTER 一致
2. **可选细粒度 chunk = 5000**：RF = 0，更细的 per-chunk 权重但元素数翻倍
3. **chunk ≤ 2000 不可用**：RF ≥ 8，权重噪声使树拓扑偏离
4. **chunk = 50000 不推荐**：元素数太少（~10），退化为 per-alignment 却仍未完美
