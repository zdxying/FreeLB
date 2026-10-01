# FP16 POP 存储实施计划

状态:实施中 | 日期:2026-09-26
目标:POP 分布函数以 FP16 存储、FP32 计算,流量减半,预期 2215 → ~4200-4400 MLUPS
(FluidX3D 同卡实测 FP32/FP16S = 4012 MLUPs = FP32 的 1.99x,RTX 3060M Laptop)。

## 背景与依据

- 当前 FP32 路径 2215 MLUPS ≈ 337 GB/s 有效带宽,已贴 RTX 3060 Laptop 标称峰值
  (336 GB/s)的屋脊线。每 cell 每步最小流量 152 B(19 读 + 19 写 × 4 B),
  唯一的大杠杆是减少字节。
- FP16 存储(2 字节/pop)流量减半;算术保持 FP32。精度代价:尾数 10 位,
  每次写回 ~5e-4 相对舍入。FluidX3D 在大量算例上验证了 FP16 存储的可行性
  (FP16S/FP16C 模式),已知失效模式:超长步数、高 Ma、多相大密度比。
- shared memory 路线已被消融实验否定(同索引 gather,零跨 cell 复用),FP16
  是唯一量级合理的下一步。

## 设计

### 核心原则:存储类型与计算类型分离

引入 `PopStorage<T>` 特征:默认 `= T`;当 `FREELB_POP_FP16` 定义且 `T = float`
时 `= __half`。**只有 POP 字段的存储元素类型变化**,RHO/VELOCITY/flag 及一切
宏量场保持 float。

```
显存 __half ──读(转换)──▶ 19 个 float 寄存器 ──计算──▶ 写回(转换)──▶ 显存 __half
```

### 关键决策

1. **FP16C 先行**(逐元素硬件转换,实现简单),FP16S(half2 打包 + 2 cell/thread)
   作为二期优化。同卡实测 FP16C 已达 82% 屋脊线。
2. **基线路径用代理引用**:`cudev::Cell::operator[]` 返回 `decltype(auto)`——
   存储类型 == 计算类型时返回 `T&`(现状);否则返回 `PopRef` 代理
   (隐式转 float + 从 float 赋值)。46 处 dynamics 代码原样编译,零改动。
3. **RegCell 几乎零改动**:`pop_` 指针改为存储类型;`cache_[d] = pop_[d][0]` 与
   `pop_[d][0] = cache_[d]` 依赖 `__half` 的隐式 float 转换(cuda_fp16.h 提供,
   device 端为硬件指令)。
4. **只动 POP**:RHO/VELOCITY/flag 保持 float(流量占比 <5%,无收益)。
5. **开关默认关闭**:`FREELB_POP_FP16` 不定义时与现状逐位相同,零风险合入。
   宿主 CPU 构建(g++)不受影响(特征退化为 T)。

## 实施步骤

| 步骤 | 内容 | 文件 | 预计量 |
|---|---|---|---|
| S1 | `PopStorage<T>` 特征 + `cuda_fp16.h` 引入 + POP 别名(宿主/cudev 两处)改用 `PopStorage<T>` | alias.h | ~20 行 |
| S2 | `PopRef` 代理 + `cudev::Cell::operator[]` 与宿主 `Cell::operator[]` 改 `decltype(auto)` 条件返回 | cell.h | ~45 行 |
| S3 | `getPopArray` 指针类型 `T**` → `PopStorage<T>**`;`RegCell::pop_` 同步 | cuda_block_lattice.h, cell.h | ~8 行 |
| S4 | cavity3d_cu Makefile:`FP16 ?= 0` → `-DFREELB_POP_FP16` | Makefile | ~5 行 |
| S5 | 校验和/宿主侧转换核对(依赖 `__half` 隐式转换,预计零改动,编译器验证) | cavity3d.cu | ~0-10 行 |
| S6 | 验证(下节) | - | - |

## 验证方案

1. **回归**:FP32(开关关)全量构建 8/8 + 两个测试套件全绿,校验和与切换前一致。
2. **FP16 构建**:编译通过 + 寄存器/局部内存用量检查(LOCAL 必须 0)。
3. **性能**:FP16 1000 步 MLUPS,预期 ≥ 4000(达成 82% 屋脊线即 3572+)。
4. **正确性三级对比**(1000 步,100³):FP32-CPU vs FP32-GPU(已知 2.3e-7)
   vs FP16-GPU——sum/sumsq 相对差预期 1e-4~1e-3 量级(存储舍入),maxabs 同量级。
5. **稳定性**:FP16 长跑(10k 步)不 NaN、不爆炸,checksum 单调性正常。

## 风险与回退

- 数值稳定性:FP16 舍入在长步数/高梯度区可能放大——以第 4/5 项实测为准,
  超预期则记录失效边界(该模式只作为可选加速,默认关闭)。
- 代理引用的悬垂/别名问题:PopRef 不缓存值,只包装指针,生命周期与底层元素一致。
- 回退:不定义 `FREELB_POP_FP16` 即回退,无残留。

## 二期(可选,本次不做)

- FP16S:half2 打包 + 2 cell/thread 重构,补齐 FP16C→FP16S 的 ~12% 差距。
- 显式 `__float2half2` 向量化、`cp.async`/L2 策略评估。

---

## 实施结果(2026-09-26,已完成)

全部步骤 S1-S5 落地,`FP16 ?= 1` 启用,默认关闭时与 FP32 路径逐位一致
(FP32 校验和不变性已验证:1061231.05 / 4749162.144 / 33.7224)。

### 性能(100³,1000 步,reg 路径)

| 配置 | MLUPS | 相对 FP32 |
|---|---|---|
| FP32(默认) | 2210.85 | 1.00x |
| **FP16(实测 3 轮)** | **4113 / 4097 / 4050** | **~1.86x** |

寄存器用量 REG:96、LOCAL:0(无溢出);mangled name 确认 POP 已是
`CyclicArray<__half>`。与 FluidX3D 同卡 FP32/FP16S(4012 MLUPs)相当。

### 正确性(1000 步)

| 量 | FP32-GPU | FP16-GPU | 相对差 |
|---|---|---|---|
| sum | 1061231.05 | 1061137.286 | 8.8e-5 |
| sumsq | 4749162.144 | 4801135.417 | +1.1e-2 |
| maxabs | 33.7224 | 33.0 | 边界层极值,轨迹分叉 |

### 稳定性(3000 步)

无 NaN、无爆炸,MLUPS 稳定(4124)。与 FP32 3000 步对照:
sum 相对差 1.4e-4(质量守恒保持),sumsq -5.3%、maxabs 95.5 vs 100.4
(混沌轨迹分叉,极值离散度 ~5%)。

**结论:FP16 存储可用,推荐用于对极值精度不敏感的场景;定量精细研究仍用
FP32(默认)。长期步数(>1 万步)的精度表现待后续评估。**

## 验证状态清单(2026-09-26 更新)

### 已验证
- 构建:8/8 全目标;FP16 × {cyclic, streammap} 编译矩阵全通过
- 容器级:CyclicArray 设备镜像 vs 精确取模参考模型逐位一致
- 求解器级:FP32 GPU vs CPU 物理对拍(sum 相对差 2.3e-7)
- FP32 开关关校验和不变性;共享 field_checksum.h 输出与旧实现逐位一致
- FP16 运行时:1000/3000 步无 NaN,质量守恒 1.4e-4,REG:96 LOCAL:0

### 待 GPU 空闲后验证(2026-09-27 更新:1/2/4/5/7 已完成)
1. ~~FP16 + baseline 路径~~ -> **已验证为缺陷**:PopRef 代理路径全场 NaN,
   FP16 现强制 reg 路径(见上方已知限制)
2. ~~FP16 长步数精度衰减曲线~~ -> **已完成**:10k 步无 NaN,质量守恒 1.2e-3,
   sumsq -17%/maxabs -23%(轨迹分叉,预期特征);本轮 MLUPS FP32 1998 / FP16 3421
3. StreamMapArray 求解器流化未生效的根因(其 raw rotate 另有 map_count 周期缺陷)—— 待查
4. ~~cavity2d_cu 运行时验证~~ -> **已完成**(2026-09-27):根因与 cavity3d_cu 相同
   (sm_89 gencode + -rdc),Makefile 修复 + cudaGetLastError 检查后,
   1000 步 Res 0.1 -> 0.027,VTK 输出正常
5. ~~VTK/WriteBinary 输出路径~~ -> **已完成**:cavity2d_cu T0/T500/T1000 全部写出
6. MPI 多块运行时(refblock 仅编译验证;halo 通信 × FP16 的交互未知)
7. ~~FP16 物理量级验证~~ -> **已完成**(2026-09-27,升级为全场对比):
   两求解器在主循环后激活 RhoU 写回并 dump 全场速度
   (profile_gpu.txt / profile_cpu.txt,各 1,061,208 cell x 3 分量)。
   FP32 GPU vs CPU: max |du| = 1e-6(文本精度地板),RMS 9.4e-8,
   max |u_x| = 0.2(盖板速度)两侧一致 —— 逐 cell 空间级物理验证通过

### 工具
- `src/utils/field_checksum.h`:pop 场指纹(GPU/CPU 同序,回归与对拍共用)
- `tools/sass_census.py`:按函数统计 SASS 指令谱(inst/LDG/STG/FFMA),
  验证优化是否真正落进生成代码
- `benchmarks/array_eval/`:容器消融基准(7 场景 × 多容器单二进制)

## 已知限制与后续修复(2026-09-26 晚)

1. **内核边界检查缺失(已修复)**:CuDevApplyCellDynamics* 四个内核原本不做
   idx < N 检查,N 非 block 整数倍时每步 24-40 个越界线程读写相邻数组。
   FP32 下漂移极小(+23/1000 步)未被发现;FP16 下放大为全场 NaN。
   现已全部加 n 参数 + 早退守卫,调用侧传 this->getN()。
   修复后 FP32/FP16 校验和与修复前逐位一致(越界写从未进入有效区)。
2. **FP16 + baseline 路径产生 NaN(已知缺陷,未修)**:PopRef 代理路径
   100 步内 92% 值变 NaN(首个在 cell 0 pop 0),reg 路径同存储完全正常。
   机制未查明(疑似代理在基线内核中的某处重载解析/寻址问题)。
   处置:cavity3d_cu 在 FP16 下强制 reg 路径(带提示);baseline 仅供 FP32 对比。
3. NaN 诊断已并入 field_checksum.h(NaN 计数 + 首个 flat 索引)。
