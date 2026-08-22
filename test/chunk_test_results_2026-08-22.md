# chunk_wtrial — chunk 大小对 RF 精度的影响（2026-08-22 修复后）

> 运行日期：2026-08-22
> 修复内容：权重相似度统计从"仅多样性位点"改为"chunk 全范围所有位点"

## 测试目的

评估 `--chunk` 参数对 chunk_wtrial 推断树拓扑精度的影响，并在修复相似度统计 bug 后重新验证。

## 数据集

| 属性 | 值 |
|------|------|
| 文件 | `example/cat/test_full.fasta` |
| 大小 | 570 MB |
| 物种 | 10 个（Felidae 猫科） |
| 参考物种 | Felis_catus |
| 置根 | Canis_lupus_familiaris |

## 运行参数

```bash
chunk_wtrial -i fasta2ref.txt --chunk <N> -t 4 --initial-round 4 --subsequent-round 2 \
  --root Canis_lupus_familiaris --log chunk_<N>.log -o chunk_<N>.tre
```

## 结果

### RF 距离汇总

| chunk | RF vs CASTER | 修复前 RF（2026-08-01） |
|-------|:-----------:|:-----------------------:|
| 1000 | **0** | 10 |
| 2000 | **0** | 8 |
| 5000 | **0** | 0 |
| 10000 | **0** | 0 |
| 50000 | **0** | 2 |

### 每 chunk 输出树

| chunk | tree |
|-------|------|
| 1000 | `(((Prionailurus,Otocolobus)100,Felis)100,Lynx)100,Acinonyx)100,Leopardus)92.6,Caracal)100,Panthera)100,Suricata),Canis);` |
| 2000 | `((((Otocolobus,Prionailurus)100,Felis)100,Lynx)100,Acinonyx)100,Leopardus)93.1,Caracal)100,Panthera)100,Suricata),Canis);` |
| 5000 | `(((Prionailurus,Otocolobus)100,Felis)100,Lynx)100,Acinonyx)100,Leopardus)93.8,Caracal)100,Panthera)100,Suricata),Canis);` |
| 10000 | `(((Prionailurus,Otocolobus)100,Felis)100,Lynx)100,Acinonyx)100,Leopardus)95,Caracal)100,Panthera)100,Suricata),Canis);` |
| 50000 | `(((Prionailurus,Otocolobus)100,Felis)100,Lynx)100,Acinonyx)100,Leopardus)96.3,Caracal)100,Panthera)100,Suricata),Canis);` |

## 结论

1. **修复后所有 chunk 值 RF=0**：相似度统计包含保守位点后，权重估算稳定，拓扑与 CASTER 完全一致
2. **bootstrap 随 chunk 增大而略增**（Leopardus/Caracal 分支 92.6→96.3），说明更大 chunk 权重更稳
3. **默认 chunk=10000 仍推荐**：RF=0，bootstrap 95，性能与精度平衡
4. 修复前 chunk≤2000 的 RF=8~10 现在全部归零，验证了 bug 修复的有效性
