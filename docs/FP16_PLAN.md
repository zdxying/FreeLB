# GPU POP 存储与策略（RegPop/DirectPop、CSE、FP16）

状态:implemented(FP16 默认关闭;RegPop 是 `benchmarks/cavity3d_cu` 默认) | 日期:2026-10-31 | 目标:D3Q19 的 FP32/FP16 × RegPop/DirectPop 四组合在 `benchmarks/cavity3d_cu` 可用并实测,FP16 只改 POP 存储且默认关闭时 FP32 逐位不变。
退出条件:FP16S(half2 打包)落地并重测、默认开启 FP16、POP 存储格式变更,或 `cavity3d_cu` 把 `--reg`/`--base` 统一为单一策略/新增第三个 POP 策略时,存储与策略决策由 `docs/ARCHITECTURE.md` §6 的 CUDA mirror 承载,本文数字由新表取代。

## 0. 结论

- FP16 存储可用,适合对极值精度不敏感的场景;定量精细研究仍用 FP32(默认)。FP16 关闭时 FP32 路径逐位不变。
- FP16 的性能收益依赖 RegPop:RegPop + FP16 是最快组合,比 FP32 + RegPop 快 87%,比 FP16 + DirectPop 快 2.9x;DirectPop 下 FP16 几乎没有收益。
- 策略是 `cudev::Cell` 的第 4 个模板参数,同一份 CSE 产物覆盖两种策略;RegPop 的 CSE 特化被选中(FP64 19 → 0)。
- `p[Q]` 地址只解析一次是性能前提;逐方向重算的回退实测见 §7。

## 1. 存储类型与计算类型分离

### 核心原则

`PopStorage<T>`(`src/utils/alias.h`)默认 `= T`;当 `FREELB_POP_FP16` 定义且 `T = float` 时 `= __half`。**只有 POP 字段的存储元素类型变化**,RHO/VELOCITY/flag 及一切宏量场保持 float。

```
显存 __half ──读(转换)──▶ 19 个 float 寄存器 ──计算──▶ 写回(转换)──▶ 显存 __half
```

- FP32 路径 2211 MLUPS ≈ 336 GB/s 有效带宽,已贴 RTX 3060 Laptop 标称峰值(336 GB/s);每 cell 每步最小流量 152 B(19 读 + 19 写 × 4 B)。FP16 以 2 B/pop 把流量减半,算术保持 FP32。
- 精度代价:尾数 10 位,每次写回 ~5e-4 相对舍入。FluidX3D 已在其 FP16S/FP16C 模式验证可行性;已知失效模式:超长步数、高 Ma、多相大密度比。同卡 FluidX3D FP16S 4012 MLUPS,本实现 FP16C 约 4113 MLUPS,相当。
- shared memory 路线已被消融实验否定(同索引 gather,零跨 cell 复用),FP16 是唯一量级合理的下一步。
- 只动 POP:RHO/VELOCITY/flag 保持 float(流量占比 <5%,无收益)。`FREELB_POP_FP16` 默认不定义;`benchmarks/cavity3d_cu/Makefile` 的 `FP16 ?= 0`,需要时 `make FP16=1`。

### PopRef

`cudev::Cell::operator[]`(`src/data_struct/cell.h`)返回 `decltype(auto)`:存储类型 == 计算类型时返回 `T&`;否则返回 `PopRef<PopStorage<T>, T>` 代理(只包装指针,隐式转 float、从 float 赋值)。dynamics 代码无需改动。

RegPop 路径几乎零改动:`PopCache<T, Q, RegPop>::v` 为 `T`(float),经 `p[d][0]` 的 `__half` 隐式转换读写存储。该转换走 `__half` 而非 `PopRef`;两者落到同一条 SASS 指令,统一只为可读性,不值得改。

`PopRef` 必须对赋值透明,否则壁面 pop 静默丢失;失败模式与证据见 §7。

## 2. 策略是 Cell 的模板参数

### D1 — 策略是追加在模板参数表末尾的类型标签,默认值在前向声明上

```cpp
// src/data_struct/pop_policy.h
struct DirectPop { static constexpr unsigned int block_size = 32; };
struct RegPop   { static constexpr unsigned int block_size = 128; };
template <typename T, typename LatSet, typename TypePack,
          typename POPPOLICY = DirectPop> class Cell;   // 默认值只此一处
// src/data_struct/cell.h
template <typename T, typename LatSet, typename TypePack, typename POPPOLICY>
class Cell : private PopCache<T, LatSet::q, POPPOLICY> { ... };
```

默认参数只能出现在第一处声明。声明放在 `pop_policy.h`,是因为 `data_struct/cuda_block_lattice.h` 的 `CellType` 等代码直接写 3 参数形式,而 `data_struct` 不得反向依赖 `lbm`;`cell.h` 经 `pop_policy.h` 即可满足,与 include 顺序无关。

CSE 特化模式必须把 `POPPOLICY` 放在**最后**,与别名元数对齐;插在中间会与别名错位,全部特化静默失配(编译通过、数值正确、只是退回循环版)。

### CSE 精确类型匹配

CSE 生成的偏特化按**精确 cell 类型**匹配:

```cpp
template <typename T, typename TypePack, typename POPPOLICY, bool WriteToField>
struct rhoUImpl<CELL<T, D3Q19<T>, TypePack, POPPOLICY>, WriteToField> { ... };
```

偏特化匹配不做基类转换:派生 cell 类命不中,会静默落回带 `for (i < LatSet::q)` 的主模板——编译通过、数值正确,只是慢。`cudev::RegCell` 只是 `Cell<T, LatSet, TypePack, RegPop>` 的别名(`cell.h` 的 `cudev::RegCell`),不是派生类。

`generated/lbm/*.ur.h` 是所有 target 共享的一份产物;把策略硬编进产物,会让使用另一种 cell 的 target 失去特化。策略因此作成 `Cell` 的模板参数,由 csegen 转发:一次发射覆盖所有策略,体积不涨,新增策略 csegen 零改动。

| 路径 | cell 类型 | CSE 特化 |
| --- | --- | --- |
| CPU `benchmarks/cavity3d` | 宿主 `Cell` | 生效(`Makefile` 的 `-D_UNROLLFOR`) |
| `examples/cavity2d_cu` | `cudev::Cell` | 生效(`cavity2d.cu` 的 `using CELL`) |
| `cavity3d_cu --base` | `cudev::Cell` | 生效(`cavity3d.cu` 的 `using CELL`) |
| `cavity3d_cu --reg`(默认) | `cudev::RegCell` = `Cell<...,RegPop>` | 生效(`cavity3d.cu` 的 `g_UseRegCell` / `using RegCELL`) |

## 3. 设计不变量 D2–D5

### D2 — 寄存器缓存放在独立 storage 类,私有继承吃 EBO

```cpp
template <typename T, unsigned Q, typename Policy> struct PopCache {};   // 空类
template <typename T, unsigned Q> struct PopCache<T, Q, RegPop> {
  PopStorage<T>* p[Q];   // q 个元素地址,构造时解析一次
  T v[Q];                // q 个值,常驻寄存器
};
```

C++17 没有 `[[no_unique_address]]`(`make.mk` 的 `-std=c++17`),空成员也占字节,必须用私有基类;否则 DirectPop 会白吃 19 个寄存器,毁掉基线路径(A1 专为此设)。

### D3 — `getNeighbor` 恒返回 `DirectPop`

邻居按值返回,调用方通常只读一个方向,给邻居分配 19 个寄存器的缓存是纯浪费;`BounceBackMovingWall` 这类一个 kernel 内取多个邻居的写法还会成倍放大。`Cell` 内的别名 `NeighborCell = Cell<T, LatSet, TypePack, DirectPop>` 明确表达这一点,不能"顺手改成继承策略"。

### D4 — 构造时解析一次地址,保留 `p[Q]` 指针数组

`RegPop` 的构造调用 `BlockLatticeBase::getPopArray` 一次填 `p`,`flush()` 把 `v` 写回;DirectPop 下两者都是 no-op,所以一个 kernel 能同时服务两种策略。

地址**不要**在 `load`/`store` 循环里逐方向重算;FP16 上这样做的实测回退与原因见 §7。

### D5 — `--reg` / `--base` 是运行时开关,策略是编译期参数

`cavity3d.cu` 的 `g_UseRegCell` 在两个已实例化的 kernel 之间选择;`block_lattice.h` 只保留两个 kernel 模板(带/不带 flag),launch 几何由 `CELLTYPE::block_size` 决定。`.cu` 的参数解析不用重写,也保住了对照实验能力。

代价:不能从 `CELLDYNAMICS` 反推 cell 类型——`tmp.h` 里有固定元数与变参两个并存的 `TaskSelector` 定义,从任务类型反推 cell 不可靠。调用点因此显式传 `CELLTYPE`,kernel 用 `static_assert(std::is_same_v<typename CELLDYNAMICS::CellType, CELLTYPE>)` 兜底,两个 `TaskSelector` 定义各提供 `using CellType = CELL`。

## 4. CPU 路径与回退

带存储/指令的实体(`PopCache`、`PopRef`、`cudev::Cell`、`RegCell`)全部在 `cell.h` 的 `#ifdef __CUDACC__` 段内;宿主 `Cell` 与 `GenericCell` 不受影响。`PopRef` 只在存储类型与计算类型不同时实例化,CPU 下 `PopStorage<T> == T`。`pop_policy.h` 只有空标签与前向声明,host 侧可见但不产生代码。生成的 host 别名 `using CELL = Cell<T, LatSet, TypePack>` 忽略第 4 参,`POPPOLICY` 退化为自由的未使用模板参数(合法),仍匹配宿主 `Cell` 的 CSE 特化。

证据:在 `examples/cavity2d` 的 CPU 翻译单元上(g++ + `-D_UNROLLFOR`)预处理,`PopCache`、`PopRef`、`RegCell`、`cudev::Cell`、`CuDevApplyCellDynamics` 均不出现;进入 CPU 编译的只有 `pop_policy.h` 的空标签、`Cell` 前向声明与两个 `block_size` 常量。

回退规则在 `make.mk` 与 `lbm.h`:`_FLB_UR_CSE` 仅在 `-D_UNROLLFOR` 且 csegen 存在时定义;否则 `lbm.h` 只编译循环版,不包含任何 `.ur.h`。`src/ur_fallback/lbm/*.ur.h` 是手写 3 参数参考,靠默认 `DirectPop` 合法,但不在 include 路径上,只供 `tools/cse` 校验。

## 5. 验证结果

D3Q19 / `CyclicArray` / RTX 3060 Laptop (sm_86) / 100³;A1–A5 为 FP32、1000 步,FP16 行以同配置 FP32 为参照,步数见项目名。汇总自 2026-09-26 起的实测。

| 项目 | 结果 |
| --- | --- |
| 构建 | 全目标构建通过;FP16 × {cyclic, streammap} 编译矩阵通过 |
| 默认关闭 | FP32 校验和 `1061231.05 / 4749162.144 / 33.7224`,与 FP32 基线逐位一致 |
| A1 DirectPop BGK | 1712 指令 / 216 FLOP / 0 FP64 / REG:44,与基线完全相同 |
| A1 DirectPop RhoU | 680 / 76 / 0 / REG:28,与基线完全相同 |
| A2 RegPop FLOP | 314 → 226(-28%) |
| A2 RegPop FP64 | 19 → 0,特化被选中 |
| A3 RegPop 寄存器 / 代码生成 | REG:96 / LOCAL:0,无溢出;mangled name 确认 FP16 时 POP 为 `CyclicArray<__half>` |
| A4 数值 | `--reg` 与 `--base` 校验和逐位相同:`sum 1061231.284 / sumsq 4749162.209 / maxabs 33.72239685 / NaN 0` |
| A5 `moment.ur.h` 体积 | 192555 B → 195378 B(+1.5%);孪生体方案会 +200KB |
| 容器 | `CyclicArray` 设备镜像 vs 精确取模参考逐位一致 |
| 内核边界 | cell-dynamics 内核带 `n` 参数与 `idx >= n` 早退;该边界检查不改变结果,FP32/FP16 checksum 逐位一致 |
| FP32 GPU vs CPU | sum 相对差 2.3e-7;全场 max 速度差 1e-6,RMS 9.4e-8,max u_x 0.2(盖板速度) |
| FP16 1000 步 | sum 相对差 8.8e-5,sumsq +1.1e-2,maxabs 33.0 vs 33.7224(边界层极值,轨迹分叉) |
| FP16 3000 步 | 无 NaN/爆炸,质量守恒 1.4e-4,sumsq -5.3%,maxabs 95.5 vs 100.4 |
| FP16 10k 步 | 无 NaN/爆炸,质量守恒 1.2e-3,sumsq -17%,maxabs -23%;MLUPS FP32 1998 / FP16 3421 |
| 测试门 | `tests/run_tests.sh` 不构建 FP16;FP16 需手动 `make FP16=1` 跑 `tests/cuda/checkpoint` |

指纹来自 `src/utils/field_checksum.h`(含 NaN 计数与首个 flat 索引);SASS 统计用 `tools/sass_census.py`。

**结论:FP16 存储可用,适合对极值精度不敏感的场景;定量精细研究仍用 FP32(默认)。**

## 6. 性能矩阵

D3Q19 / `CyclicArray` / RTX 3060 (sm_86) / 100³ lid-driven cavity / 1000 步,每格 3 遍。单位 MLUPS。

| | RegPop(寄存器驻留) | DirectPop(流式) |
| --- | --- | --- |
| **FP32** | 2179 / 2179 / 2175 | 1396 / 1338 / 1364 |
| **FP16** | 4162 / 4113 / 3960 | 1444 / 1413 / 1348 |

- 三遍均值:FP32 DirectPop 1366 / RegPop 2178;FP16 DirectPop 1402 / RegPop 4078。
- **RegPop + FP16 是最快的组合**,比 FP32 + RegPop 快 **+87%**,比 FP16 + DirectPop 快 **2.9x**。访存仍是主要瓶颈,"把 pop 缓存进寄存器"和"把 pop 减半"叠加,不互斥。
- FP32 上 CSE 的算术收益不体现:RegPop 比 DirectPop 快 1.6x,但差距完全来自访存(CSE 只减少 28% 浮点运算,访存瓶颈下不改吞吐)。
- **DirectPop 下 FP16 几乎没有收益**(1366 → 1402,+2.6%)。它每次元素访问都重新解析容器地址,数据流量没省下来。
- 因此 FP16 的性能收益依赖 RegPop;`benchmarks/cavity3d_cu` 的默认(RegPop)是正确的默认。
- 三遍均值为 FP16 + RegPop 4078 MLUPS(FP32 2211,约 1.84x);三遍范围 3960–4162。

## 7. 已知限制与开放项

### FP16

1. **长步数精度**:FP16 舍入随步数/梯度放大。10k 步仍无 NaN/爆炸,但 sumsq -17%/maxabs -23%(§5);定量精细研究用 FP32(默认),>1 万步未表征。这是可选加速的失效边界。
2. **DirectPop + FP16 仅 CSE 路径**:循环版 `src/lbm/moment.h` 的 `latset::c<LatSet>(i) * cell[i]` 没有 `Vector<int,D> * PopRef` 重载,不带 `-D_UNROLLFOR` 时 FP16 + DirectPop 编译失败;CSE 生成版(`generated/lbm/moment.ur.h`)先把 `cell[i]` 提成 `const T`,不受影响。RegPop 两条路径都能编。
3. **PopRef 必须对赋值透明**:`collision::BounceBack`(`src/lbm/collision.h`;展开版 `src/lbm/collision_unroll.h`)写作 `cell[i] = cell[iopp]`;若 `PopRef` 只有 `operator=(CT)`,重载解析会选中隐式拷贝赋值,写入落在被丢弃的临时量,壁面反弹从未到达显存。缺少写穿 `operator=(const PopRef&)` 时,FP16 在 100 步内 95% 单元变 NaN;修复后同一算例 NaN 计数 19,151,655 → 0,checksum 与 RegPop 路径在 FP16 量化噪声内一致。FP32 不实例化 `PopRef`,零影响。
4. **MPI × FP16 未验证**:refblock 仅编译验证,GPU 多块 halo 通信 × FP16 的交互未知。退出条件:GPU 多块 MPI 运行时用例以 FP16 存储构建并通过,记录质量守恒。
5. **StreamMap 容器**:求解器流化未生效,其 raw rotate 另有 `map_count` 周期缺陷。归 `docs/GPU_MULTIBLOCK_PLAN.md`;退出条件:求解器在 StreamMap 下旋转后 checksum 正确,且 `benchmarks/array_eval` 的 non-page-aligned 用例在 FP16 + StreamMap 下通过。
6. **逐方向重算地址**:去掉 `RegPop` 的 `p[Q]`、在 load/flush 里逐方向重算,FP32 持平,FP16 慢 1.75 倍(2400 vs 4000),64 位地址 load 90 → 217;载荷减半后地址解析成为主导项。不变量见 §3 D4。

### CSE 覆盖与代码生成

- **覆盖不对称**:`src/lbm/collision_unroll.h` 与 `collisionMRT_unroll.h` 是手写常驻头,不在 `make.mk` 的 `UR_CSE_BASES` 里,仍是 3 参数 `cudev::Cell<...>` → `DirectPop`。因此 `BounceBack<RegCell>` / `BounceBackMovingWall<RegCell>` 拿不到展开特化,退回循环版;这是既有的覆盖不对称,与策略参数无关,影响限于壁面格点。结果是 `moment`/`equilibrium` 有 CSE、`collision` 没有的不对称。
- 把 collision/MRT 纳入 CSE 已实测评估、暂不实施:引擎前端不支持其主模板写法(C 数组声明、`unsigned int`),也不展开 stride-2 循环;MRT 零真实使用者,BounceBack 是 0 FLOP 的纯 swap。引擎缺口清单在 `tools/cse/Makefile`。
- `collision::BGK` 自身仍是循环版,未标 `@cse`。
- `latset::w<LatSet>(k)` 仍未折叠(equilibrium/force 多处)。device 侧它经 `Fraction::operator()` 读 `__constant__`,而该函数体含 `std::cout` + `exit(-1)`(`utils/util.h`),目前只靠 nvcc 对字面量下标做常量折叠才编得过;这是目前最可能硬性搞挂 device 编译的构造。
- csegen 的 `formatConst` 仍以 `setprecision(17)` 发射裸 double 字面量;实测 nvcc 会折掉,动力学 kernel 0 条 FP64,但没有保证。
- **注释未同步**:`cell.h` 的 `cudev::RegCell` 注释块与 `cuda_block_lattice.h` 的 `BlockLatticeBase::getPopArray` 注释仍描述"逐方向重算地址",与 `PopCache::p[Q]` 相反;以代码为准。

## 8. 二期

- FP16S:half2 打包 + 2 cell/thread 重构,补齐 FP16C→FP16S 的 ~12% 差距(同卡 FP16C 约 4113 MLUPS;FluidX3D FP16S 4012 MLUPS)。
- 显式 `__float2half2` 向量化、`cp.async`/L2 策略评估。
