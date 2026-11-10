# Partial 通讯的适用范围：`NormalCommunicate()` 何时能替代 `NormalFullCommunicate()`

状态：本判据在 3D partial vs full 逐位实测落地、或 partial 路径被删除后作废；判据已定；2D 同层已实测（结论依赖主循环顺序）；3D 的 partial/full 差异按 `INTP_COMM_OVERFLOW_BUG.md` §0.0、§18 在 2D/3D 几何、单 rank 可复现，但没有 partial 逐位一致的证据；AMR 禁止（另有独立的索引越界缺陷）｜日期：2026-10-28｜目标：给定一个新算例就能确定地用 partial 还是必须用 full，并记录已否决的替代方案

> 本文只回答一个问题：partial 什么时候与 full 逐位一致。Intp/Aver 的历史与不变量在
> `docs/INTP_COMM_OVERFLOW_BUG.md`（§18 是同一问题的早先记录），本文不复述。

## 1. 一句话结论

partial 可用的条件不是"内层读不到 ghost 的漏刷方向"，而是"ghost 的未刷新分量 D̄
永不参与任何影响结果的运算"。BGK 里后者只对 D = 全 q 成立，所以 partial 不是 full 的近似，
而是一个**依赖误差传播被抵消**的技巧。

2D、同层、一层 ghost、骨架为 `Collide → Stream → Comm` 时实测逐位一致（2×2 200 步、3×3 40
步）；仓库现有的 `Stream → Collide → Comm`（`examples/cavity2d/cavity2d.cpp:252-253`、
`tests/block/normal_comm_correctness/normal_comm_correctness.cpp:110-116`）第 3–5 步分叉；
3D、AMR、ghost 参与碰撞后又被读、任何读 ghost 的 BC 或统计一律禁止（§5）。

| 情形 | partial 可用？ | 依据 |
| --- | --- | --- |
| 单块（`Comm.Comms` 为空） | 可以，两者都是空操作 | `block_lattice.hh` 的两版循环体不执行 |
| 2D 同层 `Collide → Stream → Comm` | 实测一致，仍受 §4 限制 | §3 |
| 2D 同层 `Stream → Collide → Comm` | 禁止 | §3、§5 |
| AMR / 子循环步 | 禁止 | §4、§5 第 2 条 |
| 3D | 禁止（无一致性证据） | §8 |
| ghost 被碰撞后又被读 | 禁止 | §5 第 5、6 条 |
| hidden comm | 要一致只能用 `CommMode::Full` | §6 |

## 2. 判据

两版函数结构相同，唯一差别是方向集：partial 遍历 `getCommDirection()`，full 遍历 `k < LatSet::q`
（`block_lattice.hh` 的 `normalCommunicate()` / `normalFullCommunicate()`）。方向集由
`getCommPopDir()` 填好，MPI 侧 `DirSends`/`DirRecvs` 用同一函数。记 D = 被刷新的方向集，
D̄ = 全部 q 减 D。**充要条件**：

> partial 与 full 在整条轨迹上逐位相等
> ⟺ 每次通讯后每个 ghost 格的 D̄ 都等于 full 会写进去的值
> ⟺ 所有读取 ghost D̄ 的运算（碰撞、BC、统计、hidden 打包）的输出对它不敏感，且没有别处搬它。

**为什么只有全 q 闭合**：

- D 是半空间集合：每轴取 `c(i)[axis] > 0`（或 `< 0`，符号由邻居方向决定）的方向并集去重
  （`getCommPopDir()`，`src/parallel/communicator.h`），复算得 D2Q9 面 3 / 角 5，D3Q19 面 5 / 棱 9 / 角 12。
- `opp(D) != D`（面 3 的 opp 是另外 3 个），而纯 bounce-back 的 ghost 需要 `opp(D) ⊆ D`；
  BGK 又经 rho/u 把所有 q 混进每个输出，D̄ 进一次碰撞就污染 D。

所以真子集没有一个是闭合的；partial 的正确性只能来自"D̄ 的误差还没长到被读到"。

## 3. 顺序决定量级

条件：2D D2Q9、FP64、一层 ghost、overlap=1、无 AMR、BC 为普通 bounce-back / 移动壁（在碰撞
内执行）、120²、只比 base 环（不含 ghost）、ghost 用一次 full 交换初始化。下表"强梯度"
= 初始 pop × `1 + 0.02sin(0.2x) + 0.02cos(0.13y)`（全局坐标，两侧一致）；A/B/C 顺序的数据
来自未入库探针，仓库内没有能复跑该表的测试（§8）：

| 顺序 | 2×1 | 2×2 | 3×3 |
| --- | --- | --- | --- |
| A `Stream → Collide → Comm`（仓库现状） | 首差 step 4：706 对，0.0060；step 20：1387 对，0.533 | step 5：399 对，0.133；step 20：1387 对，0.533 | 未测 |
| A + 跨缝强梯度初值 | 未单独记录 | step 3 分叉 | 未测 |
| B `Collide → Stream → Comm`（均匀初值） | 首差 step 6：2 对，0.0167；step 50：3 对，0.133 | **200 步 0 差异** | 40 步 0 差异 |
| B + 跨缝强梯度初值 | 未单独记录 | **100 步 0 差异** | 20 步 0 差异 |
| C `Collide → Comm → Stream` | 首差 step 3：706 对，与 A 同量级 | 与 A 同量级 | 未测 |

**A 为什么漏**（2×1 逐格插桩）：`ghost-in-D = 0` 每步成立——被刷新的 3 个方向在两
进程里逐位相同。漏的是 D̄：第 2 步起 715 个槽不同，第 4 步起 base 出现差异（15128 个 pop
槽里 706 个）。第一处泄漏是 base 格 `(1,2)` 从 ghost `(0,3)` 取数，而 `(0,3)` 是该块**自己的
外边界环**、**没有任何通讯写它**；`Stream()` 是整块旋转，把缝环的 D̄ 搬到块的另一侧
（`block_lattice.hh` 的 `Stream()`）。

**B 为什么在缝上看起来精确**：缝 ghost 的 9 个槽里 8–9 个已等于邻居源格，尽管只写
了 3 个——缝 ghost 与邻居 base 格是**同一个物理格**，两者自我同步，通讯只是重写已对的值。
这不是定理，见 §8。

成本（np=2、单节点、16 核、3–5 次取中位数）：

| 配置 | partial | full | payload |
| --- | --- | --- | --- |
| 2D 120² / 10×10 | 344.7 µs/step | 363.2 µs/step（+5.4%） | 31.6 → 88.3 kB（×2.80） |
| 3D 48³ / 4×4×4 | 20.6 ms/step | 28.2 ms/step（+37.3%） | 411 kB → 1.33 MB（×3.23） |
| np=1（0 条消息） | 752.0 µs/step | 770.2 µs/step（+2.4%） | 纯拆分 + 拷贝开销 |

`waitall` 12.1 → 17.0 µs（2D）、1189 → 1685 µs（3D）。2D 总差 18.5 µs，其中 waitall 增量
4.9 µs；3D 总差 7.6 ms，waitall 增量约 0.5 ms——2D 以延迟和打包为主，3D 以带宽和打包为主。

## 4. 仍然会破的地方

- **静态检查 0 违规**：2×2 / 3×3 / 2×1 / 1×2（120²、base 内缩 1）上，内层经
  `up = id - c_k` 拉的每个 `(ghost, direction)` 都落在某个通讯的 D 里。"partial 漏了内层要读的
  方向"**不是**准确描述；漏的是 §3 里 D̄ → 碰撞 → D 这条链。
- **1×N / N×1 的域边界残差**：2×1、1×2、4×1 从第 5–10 步起 2–8 个格分叉，第 20 步
  约 0.0012 → 0.017，位置在域边界格（如 loc `(3,1)`、`(58,1)`），不在缝上。full 有同源残差：
  §18.7 末段记录 `(40.5,0.5)` 处连 `NormalFullCommunicate()` 都与单块差 3 个方向。
- **AMR 的 partial 有独立缺陷**：`DirSends`/`DirRecvs` 建在
  `overlap = (level != 0) ? 2 : 1` 的 helper 帧里（`block_geometry2d.hh`），运行时却用
  `_overlap = helper.getOverlap()`（= 1），于是 `packNormalSend`（`block_lattice.hh:215`）与
  `unpackNormalRecv`（`block_lattice.hh:231`）按这些索引越界读写，打包格也落在
  `_collideFirst` 之外（先打包后碰撞）。逐几何表与 `useHelperOlap=false` 的对照见
  `docs/INTP_COMM_OVERFLOW_BUG.md` §18.8.3。
- **3D 无一致性证据**：见 §8。

## 5. 判定清单

新算例按顺序过，任意一条命中"禁止"就停：

1. **单块？** `Comm.Comms` 为空时两版循环体都不执行（`block_lattice.hh` 的
   `normalCommunicate()` / `normalFullCommunicate()`）——可用。
2. **AMR 或子循环步？** 禁止。子循环步必须用带 `(count)` 的按层过滤重载，而 hidden 及其
   partial 变体只有一版、无 count 重载（`src/data_struct/block_lattice.h` 的
   `ApplyCellDynamicsWithHiddenComm`；`NormalFullCommunicate(count)` 有按层过滤重载，
   MPI 版由同一入口驱动）；实测按层过滤的 plain 与 hidden 第 1 步即
   `ndiff=116901 maxdiff=0.0667`。子循环步里 Intp 还读粗块 halo 的 pop+rho+u（第 6 条）。
3. **主循环顺序？** `Stream` 在前 → 禁止；`Collide → Stream → Comm` 是唯一实测
   一致的骨架（§3）。
4. **2D 还是 3D？** 3D 禁止：方向集 5/9/12，partial/full 的差异可复现
   （`INTP_COMM_OVERFLOW_BUG.md` §18），但没有一条 partial 一致性证据（§8）。
5. **ghost 环是否被碰撞？** CPU 遍历 `getVoxNum()`（`ApplyCellDynamics()`，
   `block_lattice.hh`），CUDA kernel 遍历 `0..n-1`（`CuDevApplyCellDynamicsKernel`，
   `block_lattice.h`）——ghost 一律被碰撞，D̄ 必进下一轮；只有顺序 B
   被同物理格自同步抵消（§3）。
6. **是否有人读 ghost 的任一字段？** 有则禁止：`anti_pressure` 读 `getNeighbor(opp(k))` 的
   VELOCITY（`src/boundary/bounce_back_boundary.h:232`）；残差统计遍历整个数组含 ghost
   （`getToleranceRho()`/`getToleranceU()`，`block_lattice.hh`）；AMR 的 Intp 读粗块 halo 的
   pop+rho+u；CUDA 没有 cell list，只能全索引碰撞（`CuDevApplyCellDynamicsKernel`，
   `block_lattice.h`，同第 5 条）。
7. **BC 读 `getPrevious()` 还是 `getNeighbor()`？** 前者是自读该格自己的 stream 前值
   （`src/data_struct/field.h:716-726`；`block_lattice.hh:73` 的注释写成"读邻居格"是错的），
   不构成障碍；后者读 ghost，落到第 6 条。
8. **用 hidden comm 吗？** 见 §6：要一致只能用 `CommMode::Full`。

## 6. 与 hidden comm 的互斥

`Stream()` 对每个方向做一次 `CyclicArray::rotate()`，是 O(1) 的起始指针翻转
（`block_lattice.hh` 的 `Stream()`；`src/data_struct/field.h:746-755`），不搬数据——所以 hidden 的
重叠窗口只能开在两次 `Stream()` 之间的碰撞阶段内：撞邻接环 → 打包发送 → 撞其余 → 等收 →
拆包（`block_lattice.hh:1455`、`:1471`、`:1496`、`:1543`）。partial 的 D̄ 恰好在这段窗口里被
读，没有位置可躲，所以 hidden 想要一致只有 `CommMode::Full`。

一致性：hidden + Full 与 plain + Full 逐位相同——3×3 120²、200 步、全 q 全格含 ghost
环、np=1/2，maxabs 0、first_step −1；两级 AMR 20 步 np=2/4，`ndiff=0 maxdiff=0`；
`tests/feature/hidden_comm` 覆盖四种模式并对 hidden/plain 逐位比较（checks 140, failures 0）。

把 hidden 步改成 full 的成本与分解见 `benchmarks/hidden_comm/2d`（`-D_HIDDEN_PROFILE`）。

## 7. 被否决的方案

- **全 q 以外的真子集**（只补角格、只补边/面、加大方向集）：不闭包；早期未入库探针
  显示补角格不改变结果、丢单个非静止分量会改变内层解（数字不可复现，见 §8）。
- **丢 q0**：早期探针记录丢 pop 0 得 0 差异，与 §2 的判据冲突；该记录只覆盖一个配置，
  原因未解释。降级为"未解释，不得作为带宽节省依据"（探针未入库，见 §8）。
- **初值做一次 full**：D̄ 泄漏第 2 步就发生、第 4 步进入 base；§3 的 A 顺序里 t=0 强
  梯度反而更快分叉（第 3 步）。
- **只读 ghost 的旧 D̄**：阻塞点是 §5 第 6 条的四处读者——`anti_pressure`
  （`bounce_back_boundary.h:232`）、残差统计（`block_lattice.hh` 的
  `getToleranceRho()`/`getToleranceU()`）、CUDA 全索引（`block_lattice.h` 的
  `CuDevApplyCellDynamicsKernel`）、AMR Intp 的 pop+rho+u。任一处存在，D̄ 就必须正确。

## 8. 未验证

- **3D 没有一致性证据**：方向集 5/9/12 已由 `getCommPopDir()` 复算（§2）；partial/full 的
  差异按 `INTP_COMM_OVERFLOW_BUG.md` §0.0、§18 在 2D/3D 几何、单 rank 可复现，但没有
  任何 3D 下 partial 与 full 逐位一致的运行记录。
- **自同步不是定理**：§3 的 B 观察只是 2×2、5 步的测量；8–9/9 个槽相等没有证明，B 顺序的一致性
  只在大网格小步数上验证过。
- **1×N 的域边界项未定位**：只知道分叉位置，不知道由哪条 BC 与 ghost 的哪个分量触发。
- **AMR 的 `Dir*` 越界应独立立项**：与 partial/full 无关，属索引构造与运行 `overlap` 不一致。
- **hidden 路径的非 MPI 与 `_VOX_ENABLED` 边界**见 `docs/INTP_COMM_OVERFLOW_BUG.md` §18.8.4。
- **A/B/C 顺序对比没有仓库内测试**：§3 的数字来自未入库探针；仓库内
  `tests/block/normal_comm_correctness` 只跑 2D、`Stream → Collide → Comm` 的校验和，
  `tests/feature/hidden_comm` 只比较 hidden 与 plain；成本见 `benchmarks/hidden_comm/2d`（`-D_HIDDEN_PROFILE`）。

Intp/Aver 的历史、不变量与 §16/§17 的修复记录见 `docs/INTP_COMM_OVERFLOW_BUG.md`，本文不复制。
