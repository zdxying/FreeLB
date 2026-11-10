# 通信的已知缺陷（Intp 越界、MPI 收发布局、cazsAref2d 启动、同层 partial 漏刷）

状态：退出条件为 §17.4 修复、§17.7 移交 `examples/cazsAref2d` 维护、§18.8.3/§18.8.4 关闭，届时删除本文；Intp 本地越界（§1）与 MPI Intp/Aver 收发布局（§16/§17）已修，其余未修｜日期：2026-10-29｜目标：只留判据、不变量、守护映射与未修缺陷

## 0. 现在仍然成立的结论

| 结论 | 状态 | 详见 |
| --- | --- | --- |
| Intp 模板的 `+1` 行落在块外，本地路径越界写堆 | **已修**：索引钳制到块内 | §1，`tests/block/comm_index` |
| MPI Intp 收发两侧枚举的粗单元集合不同，长度不符 | **已修**：`intp::BuildStencil` 唯一枚举 | §16，`amr_comm` 的 `intp-count` |
| Aver 两侧在不同块的坐标系里取交（`recvs=10 sends=72`） | **已修**：`aver::BuildAverStencil` | §17，`cazsAref2d` OMP 构建 |
| `cazsAref2d` MPI 输出步后 `MPI_ERR_TRUNCATE` | **未修**，未定位 | §17.4 |
| 同层 partial 只刷 `getCommPopDir()` 的方向集 | **未修**：2D 量级取决于主循环顺序，3D 未测，AMR 禁止；判据与实测见 `PARTIAL_COMM_PLAN.md` | §18 |
| AMR 几何上 partial 的 `DirSends/DirRecvs` 越界（计数随几何变化） | **未修** | §18.8.3 |
| `ApplyCellDynamicsWithHiddenComm` 没有 `count` 重载，子循环 AMR 步用不了 | **未修**：对 level 过滤的 plain 路径第 1 步 `ndiff=116901` | §18.8.4 |
| 非 MPI 构建下 hidden 步不做 halo 交换 | 已知限制，不是缺陷 | §18.8.4 |

## 0.0 缺陷清单与守护

| 缺陷 | 范围 | 状态 | 守护 |
| --- | --- | --- | --- |
| Intp 索引越界（本地） | 2D/3D 几何 | **已修**（§1） | `tests/block/comm_index` |
| `cavref2d` 堆破坏 | 本地 Intp | **已修**（§1） | 同上；`amr_comm` 回归 575 项 |
| MPI Intp 收发布局不一致 | MPI 路径 | **已修**（§16） | `amr_comm` 的 `intp-count`；`comm_index` 静态检查 |
| `cazsAref2d` 未 init MPI | 示例 | **已修**（§15.1） | 无自动测试 |
| `cazsblock2d` 未 init MPI | 示例 | **已修**（§17.5） | 无自动测试 |
| `cazsAref2d` Aver 收发长度不符 | 示例 + `SharedComm` | **已修**（§17.1） | `amr_comm` 的 `aver-count` |
| `cazsAref2d` MPI 输出步后 truncate | MPI 路径 | 未修（§17.4） | 无 |
| `cazsAref2d` 细化后间歇失败 | 运行时细化路径 | **未修**；owner: `examples/cazsAref2d`（§17.7） | 无 |
| 同层 partial 通信漏刷 halo | 2D/3D，单 rank 可复现 | 未修（§18，判据见 `PARTIAL_COMM_PLAN.md`） | `tests/block/normal_comm_correctness`（默认只报告） |

## 0.1 提交索引

| 文档节 | 提交 | 改了什么 | 守护 |
| --- | --- | --- | --- |
| §1 | `c4b2a7e` | Intp 模板的 `+1` 钳制到块内（2D/3D） | `tests/block/comm_index` |
| §2 | `6723677` | `flbcomm::checkCellIndex()`：越界即报告并 `exit(1)` | 越界路径本身无测试 |
| §14 / §16 | `05c9adf` | `intp::BuildStencil` 成为跨级 Intp 的唯一枚举 | `amr_comm` 的 `intp-count` |
| §15.1 | `a438be2` | `cazsAref2d` 补 `mpi().init()` | 无自动测试 |
| §15.2 / §17 | `05c9adf` | `aver::BuildAverStencil` 统一 Aver 收发枚举 | `amr_comm` 的 `aver-count`；`cazsAref2d` OMP 构建 |
| §17.5 | `45283e0` | `cazsblock2d` 的 MPI 接线与收尾 | 无自动测试 |

表中未列出的 §17.4、§17.7、§18.8.3、§18.8.4 是未修项，没有对应提交。

## 1. Intp 本地越界（已修）

`AddtoSharedIntpComm()`（2D `block_geometry2d.hh`，3D 同构）对每个粗单元在细块上取 2×2（3D 2×2×2）邻域，模板的 `+1` 行无条件向 +y/+z 走一格；交集贴细块最后一行/列时落到块外，`popIntpCommunicate()` 按索引写 `CellPop`，不检查。

**修法/不变量**：构造侧把越界坐标钳到块内（full 含 halo），边界模板退化为取边单元；不在消费侧跳过（少写一个细单元是数值改变，且无法与"几何本就不要求这个单元"区分）。本地、MPI 收发、2D/3D 都从 `intp::BuildStencil` 出。护栏是 `tests/block/comm_index`：2D/3D 全部几何 `badSrc=0 / badDst=0`。

## 2. `flbcomm::checkCellIndex()`：报告并 `exit(1)`

`src/parallel/comm_plan.h`。越界说明几何侧的枚举错了；静默跳过会把 bug 藏进数值，且无法与"几何本就不需要这个单元"区分。11 处调用点在 `block_lattice.hh`（`unpackCrossRecv`、`averCommunicate`、`popIntpCommunicate`）与 `field_struct.hh`（`MPIAverCommunicate`、`MPIIntpCommunicate`、`AllMPIAverCommunicate`、`AllMPIIntpCommunicate` 各两套重载）；越界路径本身无自动测试。

## 3. 范围：只有 Intp 绕开几何库

同层 `Comm.Comms`、`AverComm`、`IntpComm` 三类里，同层与 Aver 两侧都走同一套 AABB 取交，边界由几何保证；只有 Intp 为凑够 4/8 个源点手写 `+1`。失效条件是细块 base 被域边界裁剪，且交集落在细块最后一行/列。探针必须用 `geo.getBlock(findBlockIndex(id))`，不能用 `geo.getBlock(id)`（按容器下标取块，block 顺序与 id 不一致时给假阳性）。

## 14. MPI 的 Intp 交换：旧枚举两侧不一致（已由 §16 取代）

旧实现接收端按"细块 self ∩ 粗块 base"（65 个粗单元）、发送端按"细块 base ∩ 粗块外扩 2 层"（66 个；`ext2` 指旧代码的 `getExtBlock(2)` 外扩块，已删）枚举，两侧 `Cells` 的元素含义不同（粗 vs 细），长度天然对不上；还有只有单侧列表的块对，发送端 `startC - 1` 会下溢成巨大的 `size_t`。试过并否决的三种对齐写法（去掉发送端偏移、发送端改用 `nbaseblock`、发送端改用接收端 `getCellIdx`）均已删除，勿重复；语义决策见 §16.2。

## 15. `cazsAref2d` 的启动阻塞（已修）

### 15.1 `LoadBalancing()` 段错误：没调 `mpi().init()`

`main()` 无参、从未调用 `mpi().init()`，`MpiManager::numTasks` 保持未初始化值（调试器里是 0），`LoadBalancing(int ProcessNum = mpi().getSize())` 据此 `BlockIndex.resize(0)`，随后 `BlockIndex[0].push_back(...)` 访问空 vector。修法是 `main(int argc, char* argv[])` + `mpi().init(&argc, &argv)`（`cazsArefblock2d.cpp:193`），与 `cavref2d`/`dambreak2d` 一致；`cazsblock2d` 同样的问题在 §17.5 一并修掉。无自动测试。

### 15.2 `AddtoSharedAverComm` 的收发长度不符

构建几何时 `setRecvSendIdx` 报 `recvs=10 sends=72 ratio=7.2`：两次 `getCellIdx` 用在不同块的坐标系（粗化虚拟块 `nCblk` 的帧 vs `self` 尺度的 AABB），18 对 10 没有 1:4 关系；病态来源是细块 22×62 与粗块 42×42 的形状差异过大。修法见 §17.1。

### 15.3 一个错误写法（勿重复）

让 `VSends` 从 `Recvs` 派生（`loc / 2` 映射到 `nCblk`）：能跑完 200 步并打印 "Calculation Complete!"，但 10 个映射全部落在 `nCblk` 之外。"能跑完"不等于"正确"，索引合法性要单独验证（`tests/block/comm_index`）。

## 16. MPI 的 Intp 交换（已修）

### 16.1 复现

双 rank 跑 `tests/block/amr_comm` 的 `intp-count`：单 rank 通过，双 rank 崩或挂（ASan 报 `unpackCrossRecv` 越界；当时的 `MPIInterpolateCommunicate` 里还有 `MPI_Barrier`，不开 ASan 时会停在那里等消息，该 barrier 已在 `4430d78` 删除）。探针显示除长度不符外，还有"只有单侧列表"的块对（如 `3→2`、`4→2` 接收为 0）。

### 16.2 语义决策：模板在粗块边界之外不允许有源点

粗源一律钳制到粗块的 full 范围（含 halo），边界处退化为取边单元。这与本地路径已接受的语义一致——单 rank 参考就是这个语义，MPI 必须与它一致。外扩方案（模板在边界处仍有完整 2×2 支撑）被否决：它要求本地路径也改回外扩枚举并跨块取数，而"钳制"已由 `comm_index` 证明可测、可解释。

### 16.3 一个枚举定义，三种布局

`src/geometry/intp_stencil.h` 的 `intp::BuildStencil<D>()`（`:104`）是"哪些粗单元参与 Intp"的唯一官方定义，本地与 MPI、2D 与 3D 都从它出：

- 支撑域 = `intersection(细 full 块, 粗 base 块)`，按粗单元计数；
- 每个粗单元按 `getIntpWeight<D>()` 的行序展开 `NumStencils = 2^D` 个 stencil（2D 4 / 3D 8），每个读钳制后的 2×2（2×2×2）粗邻域、写一个细单元；
- `SrcNum = intp::SrcNum<D>()`（`:59`）是每个 stencil 的粗源数：2D 4、3D 8；`send == recv × SrcNum` 由构造保证；
- 钳制边界是 **full 块（含 halo）**。

接入 6 处调用点（2D/3D × 本地/`InitMPIIntpComm` 发送/接收），只有三种布局：

| 调用点 | 布局 |
| --- | --- |
| `AddtoSharedIntpComm`（本地） | 每 stencil `[SrcNum 粗源, 1 细目标]` |
| `InitMPIIntpComm` 发送侧 | 每粗单元 `NumStencils × SrcNum` 粗源 |
| `InitMPIIntpComm` 接收侧 | 每 stencil 1 细目标 |

调用点：`AddtoSharedIntpComm` / `InitMPIIntpComm`（2D/3D 各一套）。老的 `getExtBlock(2)`、`startC - 1`、`getCoasenedBlock` 枚举全部删除。护栏：`unpackCrossRecv`（接收 unpack）调用 `flbcomm::checkCellIndex`；发送端索引由 `BuildStencil` 钳制生成。

### 16.4 验证

| 验证 | 结果 |
| --- | --- |
| `BuildStencil` vs 修前手写枚举（`amr_comm` 几何，逐位） | 4 块对 × 5120 项相等——本地索引不变 |
| `tests/block/comm_index`（本地 11 几何 + MPI Intp 静态检查） | 44 项全过 |
| `tests/block/amr_comm` 含 `intp-count` | 575 项全绿（该路径修复前从未跑通） |
| 双 rank 恒等式 `send == recv*SrcNum` | 4 块对全部成立（含修前为 0 的 `3→2`、`4→2`） |

限制：`amr_comm` 只有 2D；3D 的 Intp 只有 `comm_index` 的静态检查，没有数值门禁。

### 16.5 不变式（给后续维护）

1. 任何对 `intp::BuildStencil` 的改动都是**数值变更**：本地单 rank 参考与 MPI 双 rank 同时改变，`amr_comm` 无法区分"一起改对"与"一起改错"，需按 AGENTS.md §4 取 per-configuration checksum。
2. 新增跨级通信路径必须复用 `intp::BuildStencil`，不得再手写枚举；`comm_index` 的 MPI 检查按"recv == send × SrcNum + 索引在界内"审计每个块对。
3. 同类问题（Aver 的 `SharedComm` 长度不符）按同一路子修——已照此执行，见 §17。

## 17. Aver 收发统一 + 两个 CA 算例的并行（已修）

### 17.1 本地 Aver 长度不符（§15.2）的修复

`src/geometry/aver_stencil.h` 的 `aver::BuildAverStencil<D>(CBase, FBase, coarseOverlap, fineOverlap, ...)` 与 §16 同构：

- 支撑域 = `intersection(粗 FULL 块, 细 base 块)`，按粗单元计数；粗侧含朝向细块的 overlap 环（交换写的就是这个环）；
- 每个粗单元的细源 = 粗化细帧（锚点 `FBase.min - 0.5*Fvox`，格距 `2*Fvox`）中该粗单元的 2×2（2×2×2）邻域，`loc*2 + off` 后**钳制到细化网格边界**（`aver_stencil.h:158-167`）；
- `SrcNum = aver::SrcNum<D>()`（2D 4、3D 8）；三种布局（本地 `[SrcNum 细源,1 粗目标]×N`、MPI Sends、MPI Recvs）同一函数产出，`send == recv * SrcNum` 由构造保证。

接入 6 处：`AddtoSharedAverComm`（2D/3D 本地）与 `InitMPIAverComm` 收发（2D/3D），调用点 `block_geometry2d.hh:792,809,914`、`block_geometry3d.hh:898,915,1020`；旧的"粗化块 + 细化帧"两帧枚举（`getCoasenedBlock` + `getRefinedCellIdx`）删除。规整几何（`amr_comm`）上新旧本地枚举逐位相等（4 对 320 项）。`cazsAref2d` 的病态几何在新枚举下 20 个块对全部满足恒等式（`bad=0`），OMP 也能越过几何构建进入主循环——但"跑完 50000 步"是**间歇**的（§17.7）。

**不变量**：支撑域止于细 base，但邻域采样 `loc*2 + off` 被钳进细化网格 `nCmesh*2 - 1`（`aver_stencil.h:158-167`），该网格覆盖 base 加一格，所以**细源槽落到 base 之外（细 halo）是正常结果**：实测每个粗块有 2–65 个细源槽在 halo。这不是缺陷，但任何"细源一定都在 base 内"的推理都是错的。

### 17.2 overlap 不是级别规则

- `lv != 0 ? 2 : 1` 只出现在直接构造路径与 `Init(GeoHelper)`（`block_geometry2d.hh:128-131,186`）；helper 路径（本仓库全部用例）下所有块 overlap = helper 的 `_Overlap` = 1。按级别规则传参会把 `amr_comm` 的 aver/intp 搞挂（17028 越界），不要这样做。
- overlap 必须取**本块** `Block2D::getOverlap()`；Aver 与 Intp 的 MPI 分支两侧各传自己块的真实值。Intp 的 `BuildStencil` 因此是双 overlap 签名（fine/coarse 各一个）。
- 同一个"级别规则 vs 真实值"的错位在 §18.8.3 的 partial `Dir*` 索引上又出现一次，那次没有回退——它是个未修缺陷。

### 17.3 细化邻域钳制（17028 的根源，已修）

Aver 的细化邻域 `loc*2 + off` 在支撑域贴边时越过细化网格（17028 是越界索引），已钳制到 `nCmesh*2 - 1`。对内部 stencil 是空操作（逐位对比探针确认规整几何不变）。

### 17.4 `cazsAref2d` MPI 的剩余问题（未修）

修复后 MPI 双 rank 能越过旧死点并跑过输出步、Solid% 与 OMP 同量级，但在输出步之后的 field 层 Aver/Intp 交换（`BlockFieldManager::runCrossPlan`）报 `MPI_ERR_TRUNCATE`。Aver/Intp 的 20 个块对收发计数已用探针验证全部满足恒等式；truncate 来自 count 过滤（`activeDepth`）下收发双侧激活子集与请求配对的某一环，尚未定位。这是 MPI + AMR + 本算例的最后一个已知问题。

### 17.5 `cazsblock2d` 的 MPI 并行（顺带修复）

`main()` 加 `mpi().init`（同 §15.1），几何从直接构造 `BlockGeometry2D(Ni, Nj, Thread_Num, ...)`（无 rank 分配，MPI 下每 rank 块集为空 → 空跑）换成 Helper + `AdaptiveOptimization + LoadBalancing`。实测：OMP 2000 步正常；MPI np=2 完整跑完 2000 步。已知差异：块缝处 CA 捕获序列不同导致 16/91204 格的 State 差异（Conc 同位差 0.19），属 AMR 块缝上 CA 算法的串行假设，非通信错误。

### 17.6 验证汇总

| 验证 | 结果 |
| --- | --- |
| `aver::BuildAverStencil` vs 旧本地枚举（规整几何，逐位） | 4 对 320 项相等 |
| `cazsAref2d` 几何 20 块对收发恒等式（探针，np=2） | bad=0 |
| `amr_comm` aver-count / intp-count（np=2 vs 单 rank） | 逐位一致 |
| `cazsAref2d` OMP 50000 步 | 间歇跑完（§17.7）；旧代码连几何都构建不出 |
| `cazsblock2d` OMP 2000 步 / MPI np=2 2000 步 | 均完整跑完（块缝 16 格差异已记录） |

限制：3D Aver/Intp 仍只有 `comm_index` 的静态检查，无数值门禁；`cazsAref2d` 的 OMP 运行只证明能启动，数值问题见 §17.7；MPI 的 truncate 未修（§17.4）。

### 17.7 运行时细化的两个既存缺陷（未修；owner: `examples/cazsAref2d`）

owner 是 `examples/cazsAref2d` 的 AMR 运行时细化路径，与通信无关；`cazsAref2d` 因此不能当门禁：

- **16 线程，不确定**：`Thread_Num=16`、`TotalStep=50000` 下多次运行只有部分跑完，失败都在细化（第 1000 步）之后：`BlockZhuStefanescu2D::UpdateDeltaFs()` 报 `deltaf > 1`（`zhu_stefanescu2d.hh:125`）或 `Ceq > 1`，两个块的行文在 stderr 交错（`Apply_SimpleCapture` 是 `#pragma omp parallel for` over `BlockZS`，`zhu_stefanescu2d.h:423-424`），同一二进制两次运行落在不同分支——竞态。把 `RefineCheckStep` 设为 1e9（关掉动态细化、保留初始两级几何）后 3/3 全过：主循环与初始跨级通信干净，问题在"运行时细化 + 重初始化"。
- **≤8 线程，必失败且确定**：1/2/4 线程共 6 次全部在同一位置退出，`[ERROR] AABB<T, D>: min > max`，来自 `DynamicBlockLatticeHelper2D::PopFieldInit()`（`block_lattice.hh:2081`）——用**新**块号索引旧块列表：`for (iblock < BlockGeo.getBlockNum())` 之后取 `BlockGeoHelper.getAllOldBasicBlock(iblock)`（`:2104-2105`），探针实测越界那次 `n_old = 12`、`n_new = 22`，越界读到的 `BasicBlock` 再取 AABB 就 `exit(1)`。16 线程时两次细化的块数恰好相同（32/32）才绕开。

## 18. 同层级 partial 通信：只刷流入方向（未修）

`NormalCommunicate()`（`block_lattice.hh:1816`）对每个 halo 格只刷新 `getCommPopDir()`（`src/parallel/communicator.h:64`）给出的方向集 D；`normalFullCommunicate()`（`:1835`）刷全部 q。两者只差 halo 的**未刷新槽 D̄**（全部 q 减去 D）。缺陷与 AMR 无关，`np=1` 即可复现，2D 与 3D 都有；它是否可接受取决于主循环顺序、维数与 ghost 读者，判据、实测、被否决的子集与 full 的成本（机器、命令、中位数）见 `docs/PARTIAL_COMM_PLAN.md` §2–§7，这里只留缺陷本身与未修项。

受影响调用点：`cazsArefblock2d.cpp:398,406`、`openbd2d.cpp:189`、`cavity2d.cpp:253`、`cavity3d.cpp:257`（hidden 默认 Partial，输出步与末次 checkpoint 用 full：`cavity2d.cpp:266,284`、`cavity3d.cpp:270,288`）；`couette2d.cpp:215`、`cavref2d.cpp:228` 用 full，不受影响；`cazsblock2d` 的 field 路径（`cazsblock2d.cpp:424`）也不受影响。

守卫：`tests/block/normal_comm_correctness` 比较同一分解下 partial 与 full 的内层校验和；默认只报告，`--strict` 才是通过条件。AGENTS.md §3 把它记作这个已知缺口的验收条件。它只跑单进程，`tests/feature/hidden_comm` 只证明 hidden 与 plain 一致（两者错得一样），所以 **MPI 侧的 partial vs full 没有证据**。

### 18.8 AMR 与 hidden 的未修项

#### 18.8.3 AMR 几何上 partial 的 `Dir*` 索引越界（未修）

加密几何上，partial 的 `DirSends` / `DirRecvs` 索引建在 helper 的 `overlap = (level != 0) ? 2 : 1` 帧里（`block_geometry2d.hh:131`，同文件 `:186`、`:586`、`:985` 同一规则），而运行块用 `overlap = _overlap = helper.getOverlap()` = 1（同文件 `:128-131`，`useHelperOlap = true`）。缺陷在合并基线 `ff01afd` 上逐配置相同，不是本分支引入。

`tests/block/comm_index` 的 8 个几何、np=4 下实测（`N` = 块格数，计数随几何变化；以本表为准，早期未入库探针的计数不可复现，已作废）：

| 配置 | send ≥ N | recv ≥ N | 打包格在 `_collideFirst` 外 |
| --- | --- | --- | --- |
| lv2 inset | 39 | 70 | 148 |
| lv2 at corner | 6 | 6 | 55 |
| lv3 inset | 41 | 69 | 145 |
| uniform / ring refined | 0 | 0 | 0 |

不变量：越界索引由 `packNormalSend`（`block_lattice.hh:215`）读、`unpackNormalRecv`（`:231`）写；打包格不在 `_collideFirst` 里，即**先打包、后碰撞**，发出去的是碰撞前的值；同一几何上 `MPIComm` 的计数为 0，那条路径不受影响。用 `useHelperOlap = false` 重建运行几何后两侧集合一致，原因钉死在 overlap 帧错位上——与 §14（Intp 收发布局）、§15.2（Aver 收发长度）同类：两侧不在同一个帧里枚举同一个集合。

#### 18.8.4 hidden 路径的边界（未修）

- `ApplyCellDynamicsWithHiddenComm` 只有一个重载（`block_lattice.h:502`），没有 `count` / level 过滤参数，**子循环的 AMR 步不能用它**：对 level 过滤的 plain 路径实测第 1 步就 `ndiff = 116901, maxdiff = 0.0667`。
- 非 MPI 构建下 hidden 步不做任何 halo 交换，`block_lattice.hh:1555` 的 `#else` 直接退回 `ApplyCellDynamics`。写单进程对比时要知道这一点。
- `_VOX_ENABLED` 下两个半趟都会错：它们用裸 `id` 和 `flagarr[id]`（`block_lattice.hh:202-204`、`:294-296`），而 `ApplyCellDynamics` 用的是 `map[id]` / `getVoxNum()`（`block_lattice.hh:133-150`）。
- `block_lattice.hh:73` 的注释说 `getPrevious()` 读邻居格。它读的是本格 stream 前的值（`CyclicArray::getPrevious`，`field.h:716-726`），是**自读**；不要把 BC 里的 `getPrevious()` 当成"读 ghost"的来源。
