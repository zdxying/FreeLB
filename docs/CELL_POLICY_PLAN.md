# cudev::Cell 策略参数化 —— 让 CSE 生成的代码在 RegPop 路径上也生效

状态: **已实施并验证**
关联分支: dev-gpu

---

## 0. 实施结果

D3Q19 / FP32 / `CyclicArray` / RTX 3060 (sm_86),`benchmarks/cavity3d_cu` 100³ lid-driven
cavity,1000 步。基线是改动前重新编译的同配置二进制(不是历史遗留的 FP16 产物)。

| 验收项 | 基线 | 改动后 | 结论 |
| --- | --- | --- | --- |
| A1 DirectPop BGK | 1712 指令 / 216 FLOP / 0 FP64 / REG:44 | **完全相同** | 通过 |
| A1 DirectPop RhoU | 680 / 76 / 0 / REG:28 | **完全相同** | 通过 |
| A2 RegPop FLOP | 314 | **226** (-28%) | 通过 |
| A2 RegPop FP64 | 19 | **0** | 通过(特化确实被选中) |
| A3 RegPop 寄存器 | REG:96 | REG:96,无溢出(`D4` 回退后) | 通过 |
| A4 数值 | — | `--reg` 与 `--base` 校验和**逐位相同** | 通过 |
| A5 生成体积 | 192555 B | 195378 B (+1.5%) | 通过 |

关键读数:**FP64 19 → 0**。循环版主模板会引入 `LatSet::InvCs2` 之类的 double 运算,
CSE 展开后消失,这是"特化确实被选中"最可靠的证据(浮点指令总数会受其它因素干扰,
FP64 不会)。

A4 明细:`--reg` 与 `--base` 都是 `sum 1061231.284 / sumsq 4749162.209 /
maxabs 33.72239685 / NaN 0`。注意**改动前的 `--reg` 与 `--base` 并不一致**
(`1061231.05` vs `1061231.284`):旧 reg 走循环版、旧 base 走 CSE 版,浮点求和次序不同。
策略化之后两条路径共用同一套 CSE 算术,差异消失。

### 性能

完整的 2×2 矩阵(FP32/FP16 × DirectPop/RegPop)见 `docs/FP16_PLAN.md`。
要点:

| | DirectPop | RegPop |
| --- | --- | --- |
| FP32 | ~1366 | ~2178 |
| FP16 | ~1402 | **~4078** |

- **FP32 上, CSE 的算术收益不体现**:RegPop 比 DirectPop 快 1.6x,但这个差距完全来自
  访存(CSE 只减少 28% 浮点运算,访存瓶颈下不改吞吐)。这与最初的预期一致。
- **FP16 上差距才放大到 2.9x**:RegPop 4078 vs DirectPop 1402。
  访存瓶颈下,"把 pop 缓存进寄存器"和"把 pop 减半"两个手段是**叠加**的——
  FP16 配 RegPop 比 FP32 配 RegPop 还快 87%,而 FP16 配 DirectPop 几乎没有收益
  (1366 → 1402),因为它每次元素访问都重新解析容器地址,数据流量并没省下来。

### D4 是错的,已回退

最初的 D4 决定砍掉 `pop_[Q]` 指针数组、改为逐方向重算地址。**这个决定是错的。**

FP32 下它看起来是划算的:-40 个寄存器,吞吐持平(2163 vs 2179,噪声内)。
但 **FP16 下慢了 1.75 倍**:2400 vs 4000 MLUPS。

原因是地址解析的开销与数据载荷成反比。数据访存两侧完全一致
(19 次 pop 读 + 19 次写),但 64 位**地址** load 从 **90 涨到 217(+141%)**;
FP16 把载荷减半,于是地址流量成为主导项。寄存器从 96 降到 56 换来的
occupancy 收益,抵不上多出来的 127 次依赖式地址解析。

已回退为"构造时解析一次地址并保留指针数组"(即原 `RegCell` 的做法)。
回退后 FP16 的 RegPop 地址 load 从 217 降到 **4**,吞吐回到 4178 MLUPS,
寄存器回到 REG:96 / LOCAL:0(无溢出)。

**教训:局部优化的收益必须在所有配置上验证。** 只在 FP32 上量过的"省 40 个寄存器"
在 FP16 上是 1.75x 回退,而 FP16 恰恰是这个 solver 的主目标配置。

---

## 1. 问题

CSE 生成的偏特化键在**精确的 cell 类型**上：

```cpp
// generated/lbm/moment.ur.h
template <typename T, typename TypePack, bool WriteToField>
struct rhoUImpl<CELL<T, D3Q19<T>, TypePack>, WriteToField>{ ... };
```

`cudev::RegCell` 公有继承 `cudev::Cell`（`src/data_struct/cell.h:414`），但 C++ 偏特化匹配
不做基类转换, 所以 `RegCell` 匹配不上, 落回主模板 —— 主模板是带
`for (i = 0; i < LatSet::q; ++i)` 的循环版。

结果: `cavity3d_cu` 默认路径（`--reg`）拿不到 CSE, 而它恰好是唯一的默认路径。

三条路径的现状（已实测）:

| 路径 | cell 类型 | CSE 特化 | 依据 |
| --- | --- | --- | --- |
| `benchmarks/cavity3d` (CPU, g++) | 主机 `Cell` | 生效 | `cavity3d/Makefile:9` 有 `-D_UNROLLFOR` |
| `examples/cavity2d_cu` | `cudev::Cell` | 生效 | `cavity2d.cu:146` |
| `cavity3d_cu --base` | `cudev::Cell` | 生效 | `cavity3d.cu:199` |
| **`cavity3d_cu --reg` (默认)** | `cudev::RegCell` | **失效** | `cavity3d.cu:122` `g_UseRegCell = true` |

## 2. 为什么不能"直接换"

把 device pass 的 `CELL` 别名指向 `RegCell` 技术上可行, 但那是**替换**而非新增:
`cudev::Cell` 会失去特化, 于是 `cavity3d_cu --base` 和 `cavity2d_cu`（后者**根本没有**
RegCell 路径）全部退回循环版。

根因: `generated/lbm/*.ur.h` 是所有 target 共享的一份产物, 而不同 target 用不同的 cell 类。
csegen 生成时不知道下游选哪个, 把策略硬编进共享产物会让别的 target 吃剩的。

发两份特化（"孪生体"）可行, 但代价是生成体积 +200KB, 且**每加一个策略就要重发一份**。

## 3. 方案: 策略参数化

把 `RegCell` 从派生类变成 `cudev::Cell` 的一个模板参数, 并让 csegen 把该参数转发进特化。
一次发射覆盖所有策略, 体积不涨, 以后加策略 csegen 零改动。

### 3.1 SASS 证据: 两者优化的是正交瓶颈

**策略化之前**的 FP16 测量(`cavity3d_cu/cavity3d.exe`),同一 BGK 任务集合,
唯一差异是 cell 类型:

| | pop 读 LD.U16 | pop 写 ST.U16 | 其他 load | 浮点指令 | FP64 | REG |
| --- | --- | --- | --- | --- | --- | --- |
| `cudev::Cell` + CSE | 65 | 46 | 346 | **225** | **0** | 40 |
| `RegCell` 无 CSE | **19** | **19** | **67** | 314 | 19 | 96 |

- RegCell 赢在访存: 19 读 19 写, 正好等于 q, 理论最小值。
- CSE 赢在算术: 少 28% 浮点指令, 且 0 条 FP64（RTX 3060 是 GeForce, FP64 为 1/64 速率）。
- RegCell 的代价是寄存器: 96 vs 40。

**两者正交, 不可替代。** 策略化之后 RegPop 同时拿到了这两项优化
(见 §0 的验收表:FLOP 226、FP64 0、REG 96、无溢出)。

### 3.2 设计决策

**D1 — 策略是类型标签, 追加在 `Cell` 模板参数表末尾**

```cpp
// src/data_struct/pop_policy.h —— 标签放在 __CUDACC__ 之外,两侧对称
namespace cudev {
struct DirectPop { static constexpr unsigned int block_size = 32; };
struct RegPop   { static constexpr unsigned int block_size = 128; };
}
// src/lbm/equilibrium.h —— 前向声明,默认值只能给在这里
namespace cudev {
template <typename T, typename LatSet, typename TypePack,
          typename POPPOLICY = DirectPop>
class Cell;
}
// src/data_struct/cell.h —— 定义处不再重复默认值
template <typename T, typename LatSet, typename TypePack, typename POPPOLICY>
class Cell : private PopCache<T, LatSet::q, POPPOLICY> { ... };
```

追加在末尾且默认 `DirectPop` 的好处: `cavity3d.cu:199` 的三参数写法不用改;
CSE 特化模式里策略也放最后, 与别名元数对齐。默认值放在**前向声明**而非定义上,
是因为 `cell.h` 会 include 生成的 `.ur.h`(其 include `equilibrium.h`),
`equilibrium.h` 总是先被看到,而 C++ 只允许在第一个声明处给默认参数(见 §12)。

**D2 — 寄存器缓存放独立 storage 类, 私有继承吃 EBO**

C++17 无 `[[no_unique_address]]`（`make.mk:24` 是 `-std=c++17`）, 空成员占字节,
必须用私有基类:

```cpp
template <typename T, unsigned Q, typename Policy> struct PopCache {};   // 空类
template <typename T, unsigned Q> struct PopCache<T, Q, RegPop> {
  PopStorage<T>* p[Q];   // q 个元素地址,构造时解析一次
  T v[Q];                // q 个值,常驻寄存器
  __device__ void load()  { for d: v[d] = p[d][0]; }
  __device__ void store() { for d: p[d][0] = v[d]; }
};
```

不做这一步, DirectPop 会白吃 19 个寄存器, 直接毁掉基线路径。验收项 A1 专为此设。

**D3 — `getNeighbor` 恒返回 `DirectPop`**

按值返回的临时对象, 调用方通常只读 1 个方向。给邻居分配 19 个寄存器缓存是纯浪费;
`BounceBackMovingWall` 之类会用到它, 一个 kernel 内可能翻倍寄存器。
反直觉, 必须在代码注释里写明, 否则日后有人"顺手改成继承策略"会炸寄存器。

`cudev::Cell::getNeighbor` 这两个重载当前**直接调用点为零**;但宿主 `Cell::getNeighbor`
在 `freeSurface.h` / `bounce_back_boundary.h` 里有 26 处调用,所以这个 device 重载
显然是为将来把边界/CA 移植到 GPU 准备的,不是死代码,决策是面向未来的。

**D4 — 构造时解析一次地址,保留 `pop_[Q]` 指针数组**

原 `RegCell` 同时持有 `PopStorage<T>* pop_[Q]`（19 个指针）和 `T cache_[Q]`。
**保留这个指针数组**:曾尝试改成在 `load`/`store` 循环里逐方向重算地址(省 40 个
寄存器),FP32 上看着划算,但 FP16 下慢了 1.75 倍——详见 §0 的「D4 是错的,已回退」。
地址解析开销与数据载荷成反比,FP16 把载荷减半后地址流量成为主导项。

`Cell` 的构造在 RegPop 下调用 `lat->getPopArray(id, this->p)` 一次,
`flush()` 写回;DirectPop 下两者都是 no-op,这正是一个 kernel 能同时服务两种策略的原因。

**D5 — `--reg` / `--base` 保持运行时开关**

策略是编译期模板参数, 但 kernel 对两种策略各实例化一份, 在 `<<<>>>` 分支处选择。
模板数 4→2, 实例化数仍为 4（与今天持平）, 二进制体积不变, `.cu` 的参数解析不用重写,
且保住了以后做对照实验的能力。

代价: 不能用 `tmp` 的 trait 从 `CELLDYNAMICS` 反推 cell 类型 —— 两个 `TaskSelector`
定义（`tmp.h:243` 三元组版、`tmp.h:304` 变参版）的偏特化**互相歧义**
（`tmp.h:437` 注释已预警过）。改为显式传 `CELLTYPE` 实参, 并在两个 `TaskSelector`
各加 `using CellType = CELL;` 供 `static_assert` 兜住配错。

## 4. 验收标准

| # | 标准 | 基线 |
| --- | --- | --- |
| A1 | DirectPop 路径 SASS 完全无回归 | `CellKernel`: `REG:40, LOCAL:0`, 225 浮点指令, 0 FP64 |
| A2 | RegPop 路径拿到 CSE 算术 | 浮点指令降到 ~225, FP64 归零 |
| A3 | RegPop 寄存器不高于当前 | `REG:96 / LOCAL:0`，不高于改动前 |
| A4 | 数值一致 | `PopChecksumDevice` 与 CPU 基线吻合 |
| A5 | 体积不涨 | `moment.ur.h` ≈ 192KB（孪生体方案会 +200KB） |

## 5. 改动清单

### A. csegen (`third_party/cse/plugins/freelb/ur_emit.cpp`)

**A-1 preamble(约 `:190-203`)** 两侧别名元数必须一致, 否则同一份发射文本无法同时解析:

```cpp
#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack, typename POPPOLICY>
using CELL = cudev::Cell<T, LatSet, TypePack, POPPOLICY>;
#else
template <typename T, typename LatSet, typename TypePack, typename POPPOLICY>
using CELL = Cell<T, LatSet, TypePack>;   // 接受但忽略 POPPOLICY
#endif
```

host 侧忽略第 4 参后, 设备侧 `POPPOLICY` 是真策略、host 侧是自由的未使用模板参数（合法）。

**A-2 CellType shape(约 `:300-315`)** 策略放末尾, 与别名第 4 位对齐:

```cpp
// 改前
header = "template <typename T, typename TypePack" + extraDecl + ">\n";
header += "struct " + sd.name + "<CELL<T, " + latT + ", TypePack>" + extraArg + ">{\n";
aliases = "using CELLTYPE = CELL<T, " + latT + ", TypePack>;\n";

// 改后
header = "template <typename T, typename TypePack, typename POPPOLICY" + extraDecl + ">\n";
header += "struct " + sd.name + "<CELL<T, " + latT + ", TypePack, POPPOLICY>" + extraArg + ">{\n";
aliases = "using CELLTYPE = CELL<T, " + latT + ", TypePack, POPPOLICY>;\n";
```

> **全案最容易写错的一行**: 模式里 `POPPOLICY` 必须在**最后**
> (`CELL<T, latT, TypePack, POPPOLICY>`), 插在中间会与别名元数错位, 全部失配且**不报错**。

**A-3 Cell shape(约 `:295-301`)** 同样处理。当前无人使用, 但保持两条路径一致。
`TLatSet`（force.ur.h）和 `TLatSetD` 不涉及 CELL, 不动。

`equilibrium::SecondOrderImpl` 无 `using CELLTYPE`、方法不收 cell, 只有模式要改。

### B. FreeLB

- **B-1 `src/data_struct/cell.h`** — 新增 `DirectPop`/`RegPop`/`PopCache`（D2）;
  `Cell` 加第 4 参数、私有继承 `PopCache`、ctor 按策略 load、`operator[]` 按策略分派、
  新增 `flush()`（DirectPop 下 no-op）; `getNeighbor` 两重载改返回
  `Cell<T, LatSet, TypePack, DirectPop>`（D3）; 删除 `RegCell` 类, 保留同名别名:
  `template <...> using RegCell = Cell<T, LatSet, TypePack, RegPop>;`
- **B-2 `src/data_struct/block_lattice.h`** — 4 个 kernel 塌成 2 个模板,
  `flush()` 变成无条件调用; 删 `CuDevApplyCellDynamicsRegKernel` ×2。
  block size 由策略标签携带(`POPPOLICY::block_size`),`Cell` 重新导出为
  `Cell::block_size`; 曾用 `cudev::BlockSizeFor` + `Cell::PopPolicy` 别名中转,
  后一并删除——属性应长在拥有它的标签上
- **B-3 `block_lattice.h:145-166` / `block_lattice.hh:296-331`** — 4 个声明/定义变 2 个,
  签名接 `typename CELLTYPE`, 配 `static_assert(std::is_same_v<typename CELLDYNAMICS::CellType, CELLTYPE>)`
- **B-4 `src/utils/tmp.h:243, 304`** — 两个 `TaskSelector` 各加 `using CellType = CELL;`
- **B-5 `benchmarks/cavity3d_cu/cavity3d.cu`** — `:199`/`:234` 的 `using` 声明不变（别名保兼容）;
  `:276/278/289/291` 的调用补 `CELLTYPE` 实参; blockSize 不再手传
- **B-6 `examples/cavity2d_cu/cavity2d.cu`** — 补 `CELLTYPE` 实参, 策略仍 DirectPop
- **B-7 `src/data_struct/cuda_block_lattice.h:117-119`** — `getPopArray` 的注释更新。
  (曾因 D4 一度无人调用; D4 回退后它重新成为 RegPop 构造函数的唯一入口。)

### C. 文档

- `cell.h:385-412` 的 RegCell 设计说明 → 改写成策略标签说明;
  `:406-412` 那段"静默回退"警告**必须删除** —— 策略化后该陷阱不存在, 留着会误导
- `block_lattice.h:150-155` 的 "Register-resident cell dynamics" 注释
- `docs/GPU_MULTIBLOCK_PLAN.md` §8（整节引用 `RegCell` 的 `pop_[]`/`cache_[]` 布局和
  "152 B/cell" 推算）
- `docs/FP16_PLAN.md` 的 FP16+DirectPop NaN 条目 —— **注意这不是"自然消失"**:
  策略化后 `RegPop` 下 `operator[]` 返回 `T&` 确实不经过代理, 但 DirectPop 仍走
  `PopRef`, 该路径**依然坏**。实测确认后已定位并修复, 见 §11。

## 6. 分阶段实施

**阶段 1: 只开 DirectPop, 先证明无回归**（风险最低, 验收最硬）

改 A + B 全部, 但 `RegPop` 路径先不启用。跑 `cavity3d_cu --base`、`cavity2d_cu`、CPU `cavity3d`。
A1 是本阶段全部意义: `REG` 仍须 40、`LOCAL:0`、225 浮点指令、0 FP64。
任何一项漂移就停下来查 EBO。

**阶段 2: 打开 `RegPop`, 验证 CSE 落地**

核心是确认 mangled name 里 `rhoUImpl`/`SecondOrderImpl` 的实参变成
`cudev::Cell<...RegPop>`, 且 SASS 浮点指令降到 ~225、FP64 归零。
若 FP64 没归零, 说明 A-2 策略位置写错导致特化静默失配（静默失败, 不报错）。
若 `REG` 反而涨了, 说明 D4 的地址重算没被 nvcc 接受。

**阶段 3: 清理**

删 `RegCell` 类残留、更新文档、跑体积确认（A5）。

## 7. 风险

| 风险 | 触发条件 | 抓手 |
| --- | --- | --- |
| A-2 策略位置写错 → 特化静默失配 | 插在中间而非末尾 | 阶段 2 看 mangled name + FP64 是否归零。**最可能且最难发现的错**, 编译通过、数值正确、只是变慢 |
| EBO 未生效, DirectPop 白吃寄存器 | nvcc 未 elide 空基类 | A1 |
| D4 地址重算被物化成长期寄存器 | nvcc 优化取向 | A3, 不降则回退 D4 保留 `pop_[]` |
| CSE 展开顶过寄存器预算 | 19 个 pop 全活跃 + `_cse_*` 中间量 | A3 + `tools/sass_census.py` 看 `LDL`/`STL` |
| `static_assert` 误报 | 两个 `TaskSelector` 歧义 | B-4 的 `CellType` 别名绕过 trait |

## 8. 顺带修掉的既有 bug

`cell.h:370` 第二个 `getNeighbor` 重载少传了 `Lat` 参数:

```cpp
return Cell<T, LatSet, TypePack>(Id + direction * Lat->getProjection());  // 缺第 2 个实参
```

该模板从未被实例化（全仓库零调用）, 所以一直没暴露。反正要改这个函数, 一并补上。

## 9. 实施过程中发现并修掉的问题

1. **`cudev::Cell::getNeighbor(const Vector<int, LatSet::d>&)` 少传 `Lat`**
   (计划中已预判)。全仓库零调用,模板从未实例化,所以一直没暴露。

2. **`PopCache` 的 store 一开始写成值拷贝,RegPop kernel 被整条消除。**
   初版是 `auto r = raw(d); r = v[d];`,而 `raw(d)` 返回 `T&`,`auto` 推成 `T`
   ——写的是局部副本,全局内存毫无变化,nvcc 判定 kernel 无可观察副作用,
   生成的 SASS 只剩 **16 条指令**。数值上"完全正确"(因为什么都没做),
   编译无警告,只有数指令数或看 SASS 才发现。
   改为让 store 走 `(int d, T value)` 的值传递 lambda,不再依赖引用绑定。
   **教训:这类"优化"必须配一条指令数下界断言,否则静默通过。**

3. **`tools/cse/verify_moment.py` 一直在空跑。** 它的正则要求 `TypePack>` 紧邻 `>`,
   插入 `POPPOLICY` 后匹配不到任何特化,而脚本把"找到 0 个"当作通过。
   `verify_equilibrium.py` 则直接报错暴露了问题。
   两处正则改成 `TypePack(?:, POPPOLICY)?>`(必须可选,因为它同时校验
   `src/lbm/*.ur.h` 旧参考),现在 66/66 匹配,确认非空跑。

4. `verify_moment.py` / `verify_equilibrium.py` 的正则硬编码了 3 参数特化模式,
   属于"生成器改了要手工同步"的隐性耦合,已更新并加注释。

## 10. 仍然待办

- `tests/feature/fieldtemp` 编译失败(`std::get<T>` 要求 T 在 tuple 中恰好出现一次)。
  改动前后错误数相同(均 2 个),是既有缺陷,与本次无关。
- CSE 生成的浮点字面量仍是裸 double(`formatConst` 用 `setprecision(17)`),
  例如 `0.33333333333333331 * rho`。实测 nvcc 会折掉,动力学 kernel 里 0 条 FP64,
  但这是白捡的健壮性,尚未处理。
- `latset::w<LatSet>(k)` 仍未折叠(equilibrium 17 处 / force 112 处)。device 侧它经
  `Fraction::operator()` 读 `__constant__`,而那个函数体里有 `std::cout` + `exit(-1)`
  (`src/utils/util.h:96-99`),现在只靠 nvcc 对字面量下标做常量折叠才编得过。
  这是唯一能硬性搞挂 device 编译的构造,优先级最高的遗留项。
- `collision::BGK` 自身仍是循环版(`src/lbm/collision.h:54-56`),未标 `@cse`。
- 子模块 `third_party/cse` 的改动需要单独提交并更新 submodule 指针,
  否则 FreeLB 侧无法复现新的生成结果。

## 11. 后续: cudev::Cell (DirectPop) 的 FP16 支持

原先 `DirectPop + FP16` 会在 100 步内产生 95% NaN,代码里用"FP16 强制 reg"绕开,
机制记为"未查明"。现已定位并修复,`DirectPop` 完整支持 FP16。

**机制**:`cudev::Cell::operator[]` 在存储类型 != 计算类型时返回 `PopRef<__half,float>`
**prvalue**。`collision::BounceBack` 写作 `cell[i] = cell[iopp];`
(`src/lbm/collision.h:139`,`collision.ur.h` 展开后 13 处)。重载解析选中
`PopRef` **隐式拷贝赋值** `operator=(const PopRef&)`(精确匹配)而非转换用的
`operator=(CT)`(需一次用户自定义转换),写入落在被丢弃的临时量上,从未到达显存 ——
壁面反弹的 pop 一次也没被写过。

**修复**:`PopRef` 增加写穿拷贝赋值,一行:

```cpp
__device__ PopRef& operator=(const PopRef& rhs) { *ptr = *rhs.ptr; return *this; }
```

这正是"引用代理必须对赋值透明"的标准写法。四种赋值形式中只有 `cell[i] = cell[j]`
这一种此前失效,其余三种(`= float 表达式`、`= cell[j]*k`、`= cell[j]+1`)本就正常 ——
所以这个 bug 极其隐蔽:编译无警告、字段仍在演化、只是壁面不动。

修复后同一算例 NaN 计数 19,151,655 → **0**。FP32 下 `PopRef` 不实例化,零影响
(FP32 的 `--reg`/`--base` 校验和逐位不变)。

性能矩阵与结论见 `docs/FP16_PLAN.md`。要点:RegPop 在两种精度下都优于 DirectPop,
FP16 在 RegPop 上比 FP32 快 87%。**DirectPop 现在完整支持 FP16**,不再是需要绕开的路径。

## 12. CPU 路径的安全性(静态分析结论)

### 结论:CPU 编译与运行均未受影响

用 g++ + `-D_UNROLLFOR`(CPU 求解器的实际配置)预处理整个 `freelb.h`,
共 111,056 行输出,策略化引入的符号出现次数:

| 符号 | 次数 |
| --- | --- |
| `PopCache` / `block_size` | **0** |
| `PopRef` / `RegCell` / `cudev::Cell` | **0** |
| `CuDevApplyCellDynamics` | **0** |

所有策略化代码实体都在 `#ifdef __CUDACC__` 内(`cell.h` 的 `__CUDACC__` 段)。
宿主 `Cell`(`cell.h:100-191`)与 `GenericCell`(`196-257`)在第一个 diff 块之前,未触及。
CPU 求解器目录、`make.mk`、顶层 `Makefile` 全部未改动。
`PopRef` 在 CPU 上原理性不可实例化(`PopStorage<T>` 需要 `__CUDACC__`)。

CPU pass 里 CSE 特化仍正确匹配宿主 `Cell`:宿主别名
`using CELL = Cell<T, LatSet, TypePack>;` 忽略第 4 参,`POPPOLICY` 退化为自由的
未使用模板参数。

运行时旁证:FP32 的 `--reg`/`--base` 校验和 `sum 1061231.284 / sumsq 4749162.209`
改动前后逐位相同。

### `src/lbm/*.ur.h` 不在本次改动范围

曾经的 `make -C tools/cse install` 会把生成物覆盖到 `src/lbm/*.ur.h`,
那是**整份文件替换**(moment 4508 行、force 1522 行、equilibrium 538 行),
不是 POPPOLICY 的增量差异,且把"人可维护的 fallback 参考"换成了生成物。
**已回退**,`src/lbm/*.ur.h` 保持手写 3 参数版。

回退后 fallback 路径(移除 `tools/cse/csegen`)仍可编译:旧 device 别名
`cudev::Cell<T, LatSet, TypePack>` 借助 `DirectPop` 默认值仍然有效并匹配
DirectPop 特化;RegPop 在该路径下无特化、退回循环版 —— 与改动前 `RegCell`
的状况相同,不是回退。

### 前向声明的默认值必须放在 equilibrium.h

`lbm/equilibrium.h` 前向声明 `cudev::Cell` 时需要给出 `= DirectPop` 默认值。
C++ 规定默认参数只能出现在**第一个**声明上,而 `cell.h` 会 include 生成的
`.ur.h`(其 include `equilibrium.h`),所以 `equilibrium.h` 总是先被看到 ——
默认值必须在那里,不能放在 `cell.h` 的定义上。

由于这个前向声明**不在** `#ifdef __CUDACC__` 内,若默认值缺席,CPU 侧
`cudev::Cell<T, LatSet, TypePack>`(3 参数)会因参数个数不符而硬错误,而且
只在 CPU 上出错,极难排查。

修法:新增 `src/data_struct/pop_policy.h` 存放 `DirectPop` / `RegPop`
(空标签结构,放在 `__CUDACC__` 之外,host 侧可见无代价),由
`cell.h` 与 `equilibrium.h` 共同 include;并把默认参数从 `cell.h` 的定义
移到 `equilibrium.h` 的前向声明。已用一个 CPU TU 的 `static_assert` 验证:
3 参数形式解析为 `Cell<..., DirectPop>`,与 4 参数形式一致。

## 13. 相关待办

include 图的清理（`cell.h` 不使用 equilibrium 却 include 它、4 个文件靠传递 include
拿到 `equilibrium.ur.h`、其中 3 个当前构建验证不到）已单独记录在
`docs/INCLUDE_GRAPH_TODO.md`，本文件不再重复。

## 14. 遗留的不对称

`src/lbm/collision.ur.h` 与 `collisionMRT.ur.h` 是手写的、不在 `UR_CSE_BASES`
里,未重新生成,仍是 3 参数 `cudev::Cell<...>` → 绑定到 `DirectPop`。
因此 `collision::BounceBack<RegCell>` / `BounceBackMovingWall<RegCell>`
拿不到展开特化,退回循环版。这对改动前的 `RegCell` 同样如此,**不是回退**;
影响也有限(反弹壁只占壁面格点)。但确实留下了 `moment`/`equilibrium`
有 CSE、`collision` 没有的不对称,可作为后续把 collision 纳入
`UR_CSE_BASES` 的动机。
