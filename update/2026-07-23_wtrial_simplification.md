# wtrial 简化修改报告 — 2026-07-23

## 改动概览

`weighted_trial.hpp` 从 619 行缩减为 582 行（-37 行）。

## 改动一：移除 `weight_t` 类型别名，统一使用 `cnt_t`

| 位置 | 改动 |
|------|------|
| namespace 层级 | 删除 `using weight_t = double;` |
| `StepwiseColorDefaultAttributes` | `cnt_t = double`（原 `cnt_t = cnt_type`，从模板参数变为固定 double） |
| `cnt4_t` | 固定为 `double`（原为 `conditional<cnt_type, uint, ull>`） |
| 全局 | 所有 `weight_t` → `cnt_t` |

**原因**：`cnt_t` 和 `cnt4_t` 在模板参数中不可删除（支撑 DataClasses 定义），直接固定为 `double` 即可，不需要额外的类型别名。

## 改动二：模板参数从 2 个减为 1 个

```diff
- template<typename cnt_taxon_type, typename cnt_type>
+ template<typename cnt_taxon_type>
```

**原因**：`cnt_type` 原本控制 `cnt_t`（颜色组的积累计数类型），现在固定为 `double`，不需要作为模板参数。

## 改动三：DataClasses 从 5 变体减为 3 变体

```diff
- bool/uchar, uchar/uchar, bool/ushort, uchar/ushort, ushort/ushort
+ bool, uchar, ushort
```

**原因**：只有 `cnt_taxon_t` 需要类型选择，`cnt_t` 固定为 `double`。变体 5→3。

`getStepwiseColorSharedConstData()` 的 try-catch 回退链对应从 5 次减为 3 次。

## 改动四：删除 `read()` 中 `cnt_t` 的溢出检查

删除的代码块（2 个 if）：
```cpp
if (std::same_as<cnt_t, unsigned char> && nTotalSpeciesmen >= 256) { throw; }
if (std::same_as<cnt_t, unsigned short> && nTotalSpeciesmen >= 65536) { exit; }
```

**原因**：`cnt_t` 固定为 `double`，`std::same_as<cnt_t, unsigned char>` 恒为 false，这两个 if 块是死代码。

`cnt_taxon_t` 的 3 个溢出检查（`bool`、`unsigned char`、`unsigned short`）保留。

## 改动五：删除旧 `quadPosSingle(cnt_t)` 整数重载

删除了旧的 `quadPosSingle` 二重载（接受 `array<cnt_t,4>` 参数、内部使用 `cnt4_t` 类型）。保留接受 `cnt_t`（即 `double`）版本的 `quadPosSingle` 和 `quadPos`。

**原因**：旧重载是 caster_tri.hpp 的残余，wtrial 只使用 double 版本。

## 改动六：简化 pair-weight 计算

8 个单核苷酸 pair-weight 从 `pw()` 函数调用改为**直接数组读取**：

```diff
- cnt_t const a11 = pw(1,0,0), a22 = pw(2,0,0);  // pw(n1==n2) → 同核苷酸，内部直接返回 cp[col][n1]
+ cnt_t const a11 = cp[1][0], a22 = cp[2][0];    // A pairs — 直接读 colorPairWeight
```

同类：`g11=cp[1][2]`、`c11=cp[1][1]`、`t11=cp[1][3]`、`a22=cp[2][0]`、`g22=cp[2][2]`、`c22=cp[2][1]`、`t22=cp[2][3]`。

`pw()` lambda 保留，仅用于 4 个跨核苷酸计算（`r11`、`r22`、`y11`、`y22`）。

**原因**：`pw(nuc, nuc)` 内的 if 分支在 `-O3` 已被编译时消除，但直接读取更清晰地表达了"同核苷酸不需要计算"的意图。

## 验证

- 编译通过
- cat 测试数据集上 RF 距离 = 0（与 CASTER 拓扑一致）
- ref_5 多 ref 测试通过
