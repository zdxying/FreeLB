# 通信架构评审：MPI、OpenMP 与混合模式

状态：步骤 0/1/2 与 P1 已实施（`ba22f02`、`2e859f6`、`0cc16d9`、`eb63358`、`4430d78`、`19d0818`）；normal/all 与跨层 Aver/Intp 的 plan 已统一；§11 的 Order 精确化与 tag 收口已实施（`23b8cb6`、`fff29cf`、`43f14a8`、`27a109d`）｜日期：2026-11-01｜目标：记录同层通信的现状、瓶颈与剩余改进顺序；跨层 Aver/Intp 只做对照，不取代 `PARTIAL_COMM_PLAN.md` / `INTP_COMM_OVERFLOW_BUG.md`

> 退出条件：§7 的 P0/P1 项落地并由 §10 的测试覆盖后，本文压缩为 `docs/ARCHITECTURE.md` §2 的一段现状说明；在那之前，本文是通信改动的决策依据。

## 0. 范围与证据等级

本文只回答"通信架构现在是什么样、哪里会成为瓶颈、先改什么"，不重复两个已有专题文档的结论：

- partial 与 full 何时逐位一致：`docs/PARTIAL_COMM_PLAN.md`。
- Intp 越界、MPI 收发布局、`cazsAref2d` 的未修项：`docs/INTP_COMM_OVERFLOW_BUG.md`。

证据分三类，关键结论会标注：

- **[代码]** 本次逐行阅读确认，给 `file:line`。
- **[实测]** 引用仓库文档里已有的运行记录，不重新测量。
- **[推断]** 由代码结构推出、本次没有运行验证的风险。

行号是实施前的工作树（2026-10-31）；步骤 0/1/2 之后以文件与符号名为准，§9 记录已实际运行的验证。

## 1. 架构现状

同层交换的数据流：

```
BlockGeometry::InitComm()      构建 SharedComm：发送格 = 邻居 base，接收格 = 本块 ghost
        ↓
Block::Communicator            保存 Comm / AllComm / MPIComm / AllMPIComm / DirRecvs / DirSends
        ↓
BlockLattice 构造              由 Direction 推出 CommDirection / RecvDirection / SendDirection
        ↓
BlockLatticeManager::planFor() 按 peer rank 合并 link，按 (CommMode, depth) 懒构建并缓存
        ↓
runPlan() / hidden             pack → Isend/Irecv → wait → unpack
```

关键代码：

- 计划结构与合并规则：`src/parallel/comm_plan.h:52`（`activeDepth`）、`:89-100`（`CommPlan`）、`:106-130`（`buildGroups`）。
- 拓扑结构：`src/parallel/communicator.h:437-565`。
- lattice 的 plan 缓存：`src/data_struct/block_lattice.h:341-348`；构建与执行：`src/data_struct/block_lattice.hh:1046-1066`、`:1396-1434`。
- field 的 plan 缓存：`src/data_struct/field_struct.h:137-174`；构建与执行：`src/data_struct/field_struct.hh:665-786`。

设计上有两个应当保留的性质：

1. `CommPlan` 拥有 `Buffer`（`comm_plan.h:83`），发送缓冲区活到 `MPI_Waitall` 之后。
2. 两条重建路径都清空**所有深度**的计划：`BlockLatticeManager::Init()` 的两个重载（`block_lattice.hh:765-783`、`:792-806`），`BlockFieldManager` 的 `resetCommPlans()`（`field_struct.h:164-173`）。这与 AGENTS §2 的不变量一致。

## 2. MPI 模式（同层为主）

### 2.1 设计

同层 MPI 有三套消息集：

| 模式 | 刷新方向 | 刷新范围 | 典型用途 |
| --- | --- | --- | --- |
| `Partial`（`DirRecvs/DirSends`） | `getCommPopDir()`：D2Q9 面 3 / 角 5，D3Q19 面 5 / 棱 9 / 角 12 | 1 层 ghost | 省带宽；仅在很窄的条件下与 full 一致 |
| `Full`（`MPIComm`） | 全部 q | 1 层 ghost | 通用同层交换 |
| `All`（`AllMPIComm`） | 全部 q | 完整 overlap | 需要更宽 halo 的后处理 |

`Partial` 与 `Full` 的差别是方向集，`Full` 与 `All` 的差别是空间范围。命名容易让人以为 `All` 只是"更多方向"，实际它还扩大了接收格集合。

`InitMPIComm()`（`src/geometry/block_geometry2d.hh:597-759`）不重新枚举本地邻居，而是扫描 helper 的全局 geometry，按"目标 block 是否在本 rank"拆成 `DirRecvs/DirSends`；`MPIComm/AllMPIComm` 按 helper 的跨 rank 邻居和 cell 交集构建。复用几何枚举是合理的，但把 helper 与本地 block 的构造参数一致性变成了隐式前提。

### 2.2 效率

**做得好的一点**：`buildGroups()` 把同一 peer 的所有 link 合并成一条消息（`block_lattice.hh:1311-1389`）。消息数从 `O(block × neighbor)` 降到 `O(peer rank)`，这是同层路径最重要的优化。

**主要成本**：

- **pack/unpack 是单线程的 [代码]**。`runPlan()`（`block_lattice.hh:1396-1434`）在调用线程里逐个 ref 调 `packNormalSend/packFullSend`，再逐个 group 发；接收侧逐个 group `MPI_Wait` 后逐格 unpack。`packNormalSend`（`:215-226`）按方向、再按 cell 列表做标量拷贝，列面/3D 面是跨行 stride 访问。
- **先打包完所有 send，才 post 第一批 recv [代码]**。send 循环结束后才进入 recv 循环（`:1400-1417`）。因为 recv 最终在 `Waitall(send)` 之前 post，死锁风险不存在；但握手被推迟了一个完整的打包阶段。
- **plain 路径已有重叠 [代码]**。`NormalCommunicate()` 走 `normalCommunicateOverlapped()`：先 post MPI 计划，再做共享内存拷贝，最后 wait/unpack（`4430d78`）。
- **hidden 路径有重叠，但窗口被削弱 [代码]**。碰撞拆成"邻接环 → pack/send → 其余碰撞 → wait"（`:1439-1559`），结构正确；但两个碰撞半段都是串行循环（`:193-209`、`:286-301`），且 `Isend` 到 `Waitall` 之间没有 MPI 调用，能否推进网络传输取决于 MPI 实现是否有后台 progress [推断]。

**已有实测**（复自 `docs/PARTIAL_COMM_PLAN.md` §3，np=2、单节点、16 核、3–5 次取中位数）：

| 配置 | partial | full | payload |
| --- | --- | --- | --- |
| 2D 120² / 10×10 | 344.7 µs/step | 363.2 µs/step（+5.4%） | 31.6 → 88.3 kB（×2.80） |
| 3D 48³ / 4×4×4 | 20.6 ms/step | 28.2 ms/step（+37.3%） | 411 kB → 1.33 MB（×3.23） |
| np=1（0 条消息） | 752.0 µs/step | 770.2 µs/step（+2.4%） | 纯拆分 + 拷贝开销 |

同一文档记录 `waitall` 12.1 → 17.0 µs（2D）和 1189 → 1685 µs（3D）。2D 总差 18.5 µs，其中 waitall 增量 4.9 µs，约四分之一是暴露延迟，其余是打包/拆包/共享拷贝；3D 总差 7.6 ms，waitall 增量约 0.5 ms，成本主要在带宽和打包，而不是等待。2D 下 payload 涨 2.8 倍而时间只涨 5.4%，说明延迟和打包占主导；3D 下带宽和打包才明显起来。两条都指向同一优化点：先并行化 pack/unpack 与共享拷贝，再考虑继续削减 payload。

### 2.3 可扩展性

- **MPI tag 空间 [代码]**。所有交换都使用 `comm_plan.h` 的固定槽位：同层 `DirTag/FullTag/AllTag`，跨层 `AverTag/IntpTag/AllAverTag/AllIntpTag`，field 每个 manager 从 `FieldTagBase = 1 << 14` 起领取 6 个槽位。每个 Dir link 的 `Order` 是由接收端 block id + kind/index 构成的精确 `LinkKey`；`SetupMPINbrs()` 只检查 `MPI_TAG_UB` 是否覆盖预留区间（`checkTagBudget()`），block id 不再进入 tag 空间。
- **每个 rank 的全局元数据 [代码]**。`_BlockGeometry2D/3D` 已删除，helper 只保留全局 `BasicBlock` 表和 `_MPIBlockNbrs`；`InitMPIComm()` 用接收端的 `BuildSameLevelLinks()` 按需构建 `DirRecvs/DirSends`（`block_link_meta.h`），不再为每个全局 block 保存 `Communicator` cell 列表。剩余启动成本是 `_MPIBlockNbrs` 的全局邻居构建和 `SetupMPINbrs()` 的扫描。
- **每个字段类型一套消息 [代码]**。`BlockFieldManager` 为每个字段类型独立建 plan、独立 tag 区间；CA 算例里 `VELOCITY`、`CONC`、`EXCESSC` 等各自通信，消息数按字段数相乘。
- **barrier 已清理 [代码]**：带 `count` 的交换不再先做全局 barrier；lattice 的 Aver/Intp 与 field 的跨层交换也去掉了无条件 barrier（`4430d78`）。匹配只依赖 `count` 推出的 active set，`block/amr_comm` 与 `feature/hidden_comm` 覆盖这条路径。
- **MPI 线程级别 [代码]**。`mpi_manager.hh:43` 用 `MPI_Init`，等价于 `MPI_THREAD_SINGLE`。当前所有 MPI 调用都在 OpenMP 并行区之外，实际按 funneled 工作；但进程存在多个 OpenMP 线程时，严格说应使用 `MPI_Init_thread(..., MPI_THREAD_FUNNELED)` 并检查返回值 [推断：可移植性风险，未观测故障]。

### 2.4 可维护性

- **tag 命名空间 [代码]**：所有路径已收敛到 `comm_plan.h` 的固定槽位；field 的 6 槽位由每个 manager 在构造时领取，`checkTagBudget()` 只校验预留区间。跨栈并发不再依赖 block id 位宽。
- **plan 生命周期与重建规则必须两栈一致**：`postPlan/waitSends/finishRecvs/runPlan` 已抽到 `parallel/plan_runner.h`，lattice 的 `runPlan` 与 field 的 `postFieldPlan/finishFieldPlan` 共用；`initDirPlan/initFullPlan/initFieldPlan/initCrossPlan` 已收敛到 `flbcomm::buildPlan`，manager 只提供每 block 的 `CommRef` 描述。Aver/Intp 也已计划化（`43f14a8`）。buffer 归 plan、所有 depth 在 Init 时清空，是两栈共同的不变量（§1）。
- **几何与构造顺序只能有一个事实来源**：overlap 已在 `ba22f02` 统一为 helper 的 `_overlap`；全局 `_BlockGeometry2D/3D` 已删除。`BlockLattice` 构造缓存 Direction，仍要求 `InitMPIComm()` 先于 `BlockLatticeManager` 构造（`block_lattice.hh:100-113`），靠 `BlockGeometry2D(GeoHelper)` 构造函数末尾的 `InitAllMPIComm` 保证，没有类型或断言保护。`SingleBlock_OMP` 从未被构建定义，cell 级路径不可达；同层交换另有 periodic 的 `LatticeCommSet`（`periodic_boundary.h:301-321`，`:198` 注明不支持 MPI）和 hidden 调度两个入口。

## 3. OpenMP / 共享内存模式

### 3.1 设计

- 默认并行粒度是 **block**：`#pragma omp parallel for num_threads(Thread_Num)` 遍历 `BlockLats`（`block_lattice.hh:1567-1573`）。
- 同层共享交换直接按 `SharedComm::SendRecvCells` 做 `dest[id] = src[id]`，源是邻居 base、目标是本块 ghost；没有打包，也没有 MPI（`block_lattice.hh:1815-1868`）。
- `AverageCommunicate/InterpolateCommunicate` 同样是 block 级并行；Aver/Intp 的 cell 列表是 5/9 步长，`getAverage/getInterpolation` 每次重算权重（`block_lattice.hh:1871-1967`）。
- `SingleBlock_OMP` 存在时，碰撞循环改成按 cell 并行、去掉 block 级 pragma（`:138-184`），但通信函数没有对应的 cell 级实现。

### 3.2 效率与可扩展性

- 共享交换是纯内存拷贝 + 下标计算，当前用 `Array.set(id, value)` 逐格写，没有按行/面做连续拷贝；列面是 stride 访问。
- block 级并行的扩展性取决于每 rank 的 block 数。`benchmarks/hidden_comm/2d/Makefile:6`、`benchmarks/plain_comm/2d/Makefile:6` 固定 `Thread_Num=1`；示例里常见"1 block/rank"或"少量 block/rank"。后者 OpenMP 很快没有足够迭代，前者完全没有 block 级并行；`SingleBlock_OMP` 本是为这个场景准备的，但不可达。
- `schedule` 默认 static，AMR 下不同 level 的 block 计算量不同，没有 dynamic/guided，存在负载不均 [推断]。
- `Thread_Num` 是全局 `int`，所有 pragma 用 `num_threads(Thread_Num)` 强制覆盖（`src/head.h:27-38`）；想用 `OMP_NUM_THREADS` 或按阶段调整线程数必须改代码。

### 3.3 安全性与不变量

按 block 并行交换是安全的：每个线程只写自己 block 的 ghost、只读邻居 block 的 base，源/目的集合在存储上不相交 [代码]。但这个安全性依赖"发送格必须来自邻居 base、接收格必须是本块 ghost"这一不变量，代码没有断言或类型约束；`checkCellIndex` 只检查下标在 `[0, N)` 内，不检查它是不是 ghost（`comm_plan.h:139-146`）。

## 4. MPI + OpenMP 混合

- **构建层没有一等入口 [代码]**：`make mpi` 只设 `CXXC := mpic++` 和 `-DMPI_ENABLED`（`make.mk:151-153`），`make omp` 只加 `-fopenmp`（`:147-148`）。示例要么在 Makefile 里写死 `-fopenmp`（如 `cazsAref2d`、`cazsblock2d`），要么根本不开 OpenMP（如 `cavity2d`）。`make mpi omp` 会让两个 target 各自触发一次 `all`，最后一次构建的 flag 集合取决于目标顺序，不能可靠地得到两者叠加。
- **执行层是阶段式串行 [代码]**：plain 路径是 `OpenMP 共享拷贝 → 隐式 barrier → 串行 pack → Isend/Irecv → Waitall → 串行 unpack`。OpenMP 只加速第一段和碰撞，MPI 的打包/解包在单线程里，两段没有重叠。
- **hidden 路径的 OpenMP 已补齐 [代码]**：`4430d78` 把 manager 循环和两个碰撞半段改成 block 级并行，3a（同 rank 邻居的共享内存拷贝）原本就是并行的；`_HIDDEN_PROFILE` 构建仍保持串行，让计时可以归因。剩余成本主要在 pack/unpack。
- **线程级别**：如上，`MPI_Init` 而非 `MPI_Init_thread`。当前没有在 parallel region 里调 MPI，所以不会立刻出错；但这是隐式契约，后续若把打包移进 OMP 区就会踩到。
- **barrier 与 OpenMP 隐式 barrier 叠加 [代码][推断]**：每个带 count 的交换先结束一个 OMP 区域，再进 `mpi().barrier()`，再开始串行计划。AMR 一步多个交换时，同步开销会比纯 MPI 或纯 OpenMP 都高；具体量级未测。

## 5. 同层三条路径与跨层 Aver/Intp

同层是本文重点：

- **Partial**：带宽最省，但 `PARTIAL_COMM_PLAN.md` §2 的结论是"不是 full 的近似，而是依赖误差传播被抵消的技巧"；只在 2D 同层、`Collide → Stream → Comm` 且 ghost 不被读等条件下实测逐位一致。3D 没有一致性证据，AMR 明确禁止。
- **Full**：语义闭合，payload 约为 Partial 的 2.8–3.2 倍（§2.2 实测），隐藏路径要与 plain 逐位一致只能用这个模式。
- **All**：在 Full 基础上扩大空间范围，用于后处理要读更宽 halo 的场景；本地和 MPI 两侧都维护 `AllComm/AllMPIComm`。

**默认值 [代码]**：`ApplyCellDynamicsWithHiddenComm` 的模板默认已在 `ba22f02` 改为 `CommMode::Full`；旧默认 `Partial` 会把调用者推向未验证组合（条件见 `PARTIAL_COMM_PLAN.md`）。非 MPI 构建会回落到普通 `ApplyCellDynamics`（`block_lattice.hh` 的 hidden 入口），所以这只影响 MPI 构建。

**跨层 Aver/Intp 已计划化 [代码]**：`43f14a8` 把 lattice 和 field 的 Aver/Intp 改成按 peer 合并的 plan：固定 tag、offset pack/unpack、plan 持有 buffer。发送/接收两侧仍各自做深度过滤，规则与旧的 `count % 2^deLevel` 一致；`block/amr_comm` 的 575 项检查覆盖了这些边界。

已知未修项不在本文重新立项：同层 partial 只刷流入方向与 AMR 几何上 partial 的 `Dir*` 越界（`docs/INTP_COMM_OVERFLOW_BUG.md` §18、§18.8.3），以及 `cazsAref2d` 输出步后的 `MPI_ERR_TRUNCATE` 和运行时细化的间歇失败（同文 §17.4、§17.7）。状态与守护以该文 §0.0 的清单为准。

## 6. 架构判断

1. **同层 MPI 的主结构不需要推翻**：每 peer 一条消息、plan 缓存、buffer 归 plan，这三点是正确的，应保留。
2. **下一步的收益主要在"搬数据"而不是"发消息"**：pack/unpack、共享拷贝、unpack 的串行化是实测成本的大头。
3. **OpenMP 设计与 MPI 设计没有对齐**：block 级并行在"多 block/rank"时够用，在"少 block/rank"时失效；cell 级路径不可达；通信本身没有并行化。
4. **tag 空间与全局 helper geometry 已在步骤 1 处理**：Dir/Full/All 用常量 tag + receiver-derived Order，helper 只保留 BasicBlock/邻居表；剩余可扩展性风险是 field tag 分离、每字段一套消息和 barrier。
5. **可维护性风险集中在剩余的重复**：plan post/wait/unpack、normal/all builder 与跨层 Aver/Intp 已共用；剩余是 periodic/hidden 仍是独立入口。

## 7. 建议与优先级

### P0：正确性与默认值

实施状态：P0-1/P0-3 已在 `ba22f02`/`2e859f6` 落地；P0-2 的 tag 收口在 `fff29cf`/`43f14a8`/`27a109d` 完成；P0-4 只对 hidden 路径硬拒绝（`ba22f02`），plain `NormalCommunicate(count)` 保留为 `amr_comm` 的诊断路径。

| 项 | 改动点 | 验收 |
| --- | --- | --- |
| P0-1 hidden 默认模式 | 把 `ApplyCellDynamicsWithHiddenComm` 的默认 `Mode` 改为 `Full`，或去掉默认强制显式传入；`cavity2d` 显式写出 Mode | `feature/hidden_comm` 四种模式仍通过；`cavity2d` 的 checksum 变化需解释 |
| P0-2 tag 空间 | 同层已是一 peer 一条消息，group tag 可退化为"每种交换一个小常量"，`Order` 继续做打包顺序；建 plan 时查 `MPI_TAG_UB` 并对超限报错；field Aver/Intp 也纳入同一策略 | `block/plain_comm_mpi`、`block/refblock`、`block/amr_comm`、`feature/hidden_comm` |
| P0-3 overlap 单一来源 | `Init(GeoHelper)`、两个构造函数、`InitMPIComm/InitMPIAverComm` 用同一个函数决定 block overlap | `block/comm_index`、`block/amr_comm` 的 `intp-count/aver-count` |

### P1：效率

实施状态：1/2/3/5/6 已实施；4 已由 `43f14a8` 完成（按 peer 合并 + plan 化）。

测量（2026-10-31，i7-11800H / 16 线程，OpenMPI，`OMP_NUM_THREADS=4 mpirun --bind-to none -np 4`，3 次中位数）：

- hidden 2D，120² / 10×10，300 步，full，构建不带 `_HIDDEN_PROFILE`：P1 前 `eb63358` 12.31 MLUPS（10.83–13.37），P1 后 21.34 MLUPS（19.51–21.68），约 +73%。收益主要来自 hidden 两半碰撞的 block 级并行，pack/unpack 也参与了，但未拆分归因。
- plain 2D，同配置、500 步：partial 24.38 → 24.44，full 20.69 → 19.91，all 22.32 → 20.10 MLUPS。没有稳定收益，full/all 略慢；小 plan 的并行/等待开销抵消了 recv-first 与共享拷贝重叠。不能声称 plain 提升。

1. **recv 先 post，再 pack**：`runPlan` 和 hidden 都改；先 post 全部 `Irecv`，再打包发送。
2. **pack/unpack 按 group 并行**：group 之间 buffer 不相交，可做 OpenMP；group 内保持串行。
3. **plain 路径做共享拷贝与 MPI 的重叠**：post recv → 共享拷贝 → wait；hidden 路径把两个碰撞半段改成 block 级并行。
4. **Aver/Intp 复用 plan/buffer**：至少去掉每次调用的三层 vector 分配；按 peer 合并是进一步目标。
5. **审查 barrier**：带 count 的交换不需要全局 barrier 来保证匹配；删除或缩小到 I/O/检查点前。保留时必须写清楚它匹配的是什么。删除前用 `block/amr_comm`、`feature/hidden_comm` 验证。
6. **增加 `make hybrid`**：`mpic++ + -DMPI_ENABLED + -fopenmp`；`MPI_Init` 换成 `MPI_Init_thread(MPI_THREAD_FUNNELED)` 并检查返回值。

### P2：可维护性

1. **合并 lattice 与 field 的 plan 引擎**：normal/all 与跨层 Aver/Intp 已实施——`flbcomm::buildPlan` 负责 depth 过滤、分组、offset 和 buffer sizing，manager 只提供每 block 的 `CommRef`（含收/发侧的 tag/order/size 策略）。验收：`make test-cpu`，`block/plain_comm_mpi`、`block/amr_comm`、`feature/hidden_comm` 通过，`regression/cavity2d` 无 checksum 变化。
2. **抽出唯一的 active-set 判定**：`isActive(level, count, role)`，替换多处 `std::pow(2, deLevel)` 和 `deLevel±1` 的手写过滤。验收：`block/amr_comm` 的 `aver-count/intp-count` 不变，`feature/hidden_comm` 四种模式仍通过。
3. **整理通信结构代际**：`BasicCommSet/LatticeCommSet` 仍被 periodic boundary 使用，但和 `SharedComm/DistributedComm/CommPlan` 混在同一层头文件；要么并入，要么明确标注 legacy/local-only。验收：`tools/check_headers.sh`、`feature/instantiation_smoke`、`block/refblock` 仍通过。
4. **文档同步**：`docs/ARCHITECTURE.md` §2 补一句"field 也有自己的 plan 栈与 tag 区间"；`PARTIAL_COMM_PLAN.md` 在默认值改变后更新 §6。验收：命令与链接可解析，表述与代码一致（AGENTS §3 `docs/` 行）。

## 8. 被否决或暂不采纳的方案

| 方案 | 结论 | 理由 |
| --- | --- | --- |
| 立即把两套通信栈重写成一套 | 暂不采纳 | 改动面覆盖 geometry、lattice、field、periodic 四处；先做 P0/P1 能在不动公共 API 的前提下拿到主要收益 |
| 删除 Partial，只保留 Full | 否决 | 2D 同层在 `Collide → Stream → Comm` 下有实测逐位一致和带宽收益；保留为显式 opt-in，而不是默认 |
| 用 `MPI_Neighbor_alltoallw` 取代逐 peer `Isend/Irecv` | 暂不采纳 | 每个 link 的 payload 随 level/面/角变化，还要维护本地/远程混合链路和双方一致的打包顺序；当前消息数 `O(peers)` 不是瓶颈，先解决 pack/unpack |
| `MPI_THREAD_MULTIPLE` + 并行区里直接调 MPI | 否决 | 当前所有 MPI 调用都在并行区外；FUNNELED 足够，MULTIPLE 会引入实现相关的锁开销和匹配风险 |
| 把所有 `mpi().barrier()` 一次删光 | 暂不采纳 | 匹配语义上多数可删，但 rank 间 level 过滤和 I/O 顺序需要测试确认；按 §7 P1-5 分批做 |
| 把 field 的消息并入 lattice 的同一条消息 | 暂不采纳 | 字段由不同 manager 拥有，生命周期和 pack 逻辑不同；作为后续优化记录，不阻塞当前 P0/P1 |
| 把 AMR + Partial 列入 P0（运行时拒绝或修到有测试） | 部分采纳 | hidden 路径已在 `ba22f02` 硬拒绝；plain `NormalCommunicate(count)` 保留给 `amr_comm` 诊断，根因仍是 `INTP_COMM_OVERFLOW_BUG.md` §18.8.3 的 `Dir*` 越界，修复属该文 |

## 9. 未验证与风险

- 步骤 0/1 已由 `make test-cpu` 验证（PASS 17, FAIL 0, SKIP 2）；§2.2 的性能数字仍来自 `PARTIAL_COMM_PLAN.md` 的既有实测，本次没有新的 MLUPS。
- tag 超限没有在真实大分解上复现；`checkTagBudget()` 在 geometry 初始化时校验 `[FieldTagBase, FieldTagBase + FieldTagSpace)` 能放进 `MPI_TAG_UB`，block id 已不参与 tag。
- barrier 的实际成本没有测；"可以去掉"是基于匹配语义的判断，删除前需要测试。
- hidden 碰撞串行对混合构建的实际影响没有测；代码层面可以确认没有 OMP pragma，损失取决于 block 数和 `Thread_Num`。
- hidden 的重叠是否真的推进网络传输取决于 MPI 实现；实测 `waitall` 仍有 12–17 µs（2D），说明重叠不彻底。
- 3D partial 仍无一致性证据；`FixedPeriodicBoundaryManager` 明确不支持 MPI（`src/boundary/periodic_boundary.h:198`），periodic 路径目前只有本地交换实现。
- `SingleBlock_OMP` 的"不可达"来自全仓库 grep；如果用户在命令行手动 `-DSingleBlock_OMP`，该路径会编译，但通信函数仍没有 cell 级并行。

## 10. 复现与验收

改动落地时按 AGENTS §3 的行执行；与本主题直接相关的目标：

```bash
make test-cpu
tests/run_tests.sh --no-cuda
```

通信专项（在 `tests/` 下由 `run_tests.sh` 统一构建/运行，也可单独进目录 `make`）：

- `block/plain_comm_mpi`、`block/refblock`：同层 MPI 计划。
- `block/amr_comm`：跨级 Aver/Intp 的 `aver-count`、`intp-count`。
- `block/comm_index`：cell 索引与 stencil 范围。
- `block/normal_comm_correctness`：partial/full 对比；`--strict` 目前按已知缺陷为红，记录计数。
- `feature/hidden_comm`：hidden 四种模式与 plain 的逐位比较。
- `block/checkpoint_decomp`、`feature/checkpoint`：计划重建后仍能重启。

性能复测用 `benchmarks/hidden_comm/2d`（`-D_HIDDEN_PROFILE`）与 `benchmarks/plain_comm/2d`，记录机器、命令和 MLUPS；neutral 也是有效结果。

## 11. Order 精确化与 tag 收口（2026-11-01 批准）

状态：已实施（`23b8cb6`、`fff29cf`、`43f14a8`、`27a109d`）｜目标：删除 Dir 排序键的哈希碰撞面，并让 block id 只出现在排序键、不再出现在 MPI tag 里。

> 退出条件：§11.3 的四个提交全部落地、`make test-cpu` 全绿且 `regression/cavity2d` checksum 不变后，本节压缩为 `docs/ARCHITECTURE.md` §2 的一段现状说明。

### 11.1 Order 的表示：4 字段 LinkKey

- 决定：`Order` 改为 `flbcomm::LinkKey{Receiver, Sender, Kind, Index}`，按字段字典序比较；新增 `src/parallel/link_key.h`。
- 理由：排序只在建 plan 时做一次，比较成本无关紧要；精确全序删除了“两个 link 哈希碰撞后 `std::sort` 在两端排出不同次序、静默交换 halo”这一类失败模式。排序键仍由两端从几何独立推出。
- 否决：
  - 保留 splitmix64 hash：收益只是免去位宽假设，而该假设已由 tag 检查强制；代价是一个静默错误面。
  - 64 位精确打包（receiver/sender 各 14 bit）：改动更小，但把排序键和 `FieldTagBase` 的位宽约束绑定；4 字段不依赖任何上限。
- 数值影响：排序只决定消息内 link 的字节布局，不改变格子到格子的映射；`regression/cavity2d` 的 checksum 必须保持不变，变了就是映射写错。

### 11.2 tag 策略：常量 + per-manager 基址

tag 只需在同一 communicator、同一对 rank、同时在途的消息之间唯一。同层交换已经是“每 peer 一条消息 + 常量 tag”；本节把同一规则推广到 field normal/all 与跨层 Aver/Intp。

- 固定槽位：`DirTag/FullTag/AllTag/AverTag/IntpTag/AllAverTag/AllIntpTag`。
- field plan：每个 `BlockFieldManager` 构造时从 `allocateFieldTagBase()` 取一个基址，normal/all 用基址 +0/+1；不再用 `FieldTagBase + blockId`。
- 跨层 Aver/Intp：按 peer 合并成一条消息，复用 `buildPlan`/`plan_runner`，tag 用固定槽位。
- 前提：所有 rank 以相同顺序构造相同的 field manager 集合。违反时 tag 不匹配会挂起或报错，不会静默错数据。
- 否决：保留 block id tag（block id 位数与 tag 空间耦合，且两个 field manager 并发时会撞 tag）；为每个子系统 `MPI_Comm_dup`（额外 communicator 生命周期，收益不抵成本）。
- 结果：`checkTagBudget()` 退化为只校验 `MPI_TAG_UB` 能容纳预留区间，block id 不再参与 tag。

### 11.3 改动分解

| 提交 | 内容 | 主要文件 |
| --- | --- | --- |
| 1 | `LinkKey` 替换 hash；Full/All/Field 用 `{receiver, sender, 0, 0}` | `parallel/link_key.h`（新增）、`comm_plan.h`、`communicator.h`、`block_link_meta.h`、`block_lattice.hh`、`field_struct.hh` |
| 2 | field plan 改用 per-manager tag 基址 | `comm_plan.h`、`field_struct.h/hh` |
| 3 | Aver/Intp 计划化：按 peer 合并、常量 tag、offset pack/unpack | `comm_plan.h`、`block_lattice.h/hh`、`field_struct.h/hh` |
| 4 | 收缩 `checkTagBudget()`，删除 block-id 检查 | `comm_plan.h`、`block_geometry2d.hh`、`block_geometry3d.hh` |

### 11.4 验收

- 每个提交跑 `make test-cpu`；涉及 `src/parallel` 的提交另跑 `block/plain_comm_mpi`、`block/refblock`、`block/amr_comm`、`block/comm_index`、`block/checkpoint_decomp`、`feature/hidden_comm`、`block/normal_comm_correctness`。
- `regression/cavity2d` checksum 必须不变（排序与 tag 改动不改变数值）。
- `tools/check_headers.sh` 退出 0；新增头文件纳入检查。
- Phase 3 用 `block/amr_comm` 的 `aver-count`/`intp-count` 覆盖过滤条件边界。
