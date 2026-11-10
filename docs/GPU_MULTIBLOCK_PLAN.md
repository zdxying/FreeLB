# GPU 多 block 支持评估计划

状态:评估(未实施) | 日期:2026-09-26 | 机器:RTX 3060 Laptop / sm_86 / 30 SM / 6 GB / CUDA 12.5(WSL2)

## 0. 摘要

- **现状**:GPU 路径只跑 `BlockLats[0]`,`block_num > 1` 时静默只算 1/N 的格子 —— 这是正确性隐患,不是性能问题,应该先堵。
- **必要性**:只有 **AMR(跨分辨率)** 和 **多 GPU** 这两件事真正必须依赖多 block。显存、稀疏几何这两条常见理由经核算**不成立**(见 §2),不要拿它们当立项依据。
- **性能损耗**:主导项是 kernel 启动,不是带宽也不是波前量化。本卡 WSL 实测边际 ~6.6 µs/launch;每 block 每步 2 次 launch,故 8 block 已约 -30%、16 block 约 -100%。**CUDA Graph 可把 N≤32 的损耗压回 ≤15%**,这是决定"要不要做"的关键变量。
- **建议**:先做 M0 安全护栏(零风险),再做 M1 单卡同分辨率多 block + Graph;AMR/跨分辨率另立一期,先做可行性判断再谈实现。

---

## 1. 现状勘察(代码事实)

| 位置 | 事实 |
|---|---|
| `src/data_struct/block_lattice.hh:1024` | 注释 `for now only one block is supported`;`CuDevStream` / `CuDevApplyCellDynamics` / `...Reg` 四处硬编码 `BlockLats[0]`,多 block 循环已写好但被注释 |
| `README.md:31` | `Cuda for multi-block structure is not supported for now` |
| `benchmarks/cavity3d_cu/cavity3d.ini` | `block_num = 1`(GPU 基准只测单块) |
| `src/data_struct/block_lattice.hh:292` | `CuDevStream` 每 block 一次 `<<<LatSet::q, 1>>>` launch(19 个 CTA × 1 线程),只做指针 rotate |
| `src/data_struct/cuda_field.h` (`cudev::CyclicArray`) | 块内 streaming 是**指针旋转**,O(1) 与格子数无关;跨 block 与跨分辨率不适用 |
| `src/data_struct/block_lattice.hh` (`NormalCommunicate`, `MPI*Communicate`) | 块间通信**全部走宿主内存**;无 device 端版本 |
| `src/data_struct/block_lattice.h:97` | 每 block 一套 `dev_Delta_Index` + `dev_Fields` + `devOmega/dev_Omega/dev_fOmega`(3 次 1 元素 `cudaMalloc`)+ `dev_BlockLat` |
| `src/data_struct/block_lattice.h` (`CuDevApplyCellDynamicsKernel`) | 已补 `idx < N` 边界检查;`Reg` 变体已并入同一个 kernel,POP 策略由 `cudev::Cell` 的 `POPPOLICY` 模板参数决定 |
| `src/utils/field_checksum.h` | `PopChecksumDevice` 单 block;调用处 `getBlockLat(0)` |
| `src/data_struct/block_lattice.h:64` | Omega / _Omega / fOmega 已是**每 block 独立** —— 不同分辨率用不同 τ 的结构已经具备 |
| `src/geometry/block_geometry3d.h:172` | `getMaxLevel()` 存在 → AMR 多分辨率 block 是一等公民(CPU 侧) |

**顺带发现的两个既有缺陷(与多 block 强相关,建议 M0 一并修)**

1. **越界写**:`blockNum = ceil(N/128)`,最后一个 CTA 有 `N mod 128` 个越界线程,而 Reg/non-Reg 两个 kernel 都不做边界检查就 `cell.flush()` 回写。单块时最多越界 127 个 cell(100³ 下是 64 个),**多 block 后越界次数 × block 数**。
2. **静默错误**:`block_num > 1` 时 GPU 只算 block 0,其余 block 的 POP 停在初值,程序照常打印 `Calculation Complete!` 和 MLUPS。这与 Makefile 里记录的那次 sm_89 事故是同一类问题(launch 失败不被检查)。

---

## 2. 必要性分级

### 2.1 成立 —— 没有多 block 就做不到

| 需求 | 为什么必须多 block |
|---|---|
| **AMR / 跨分辨率** | 不同加密层级 = 不同 Δx、不同 τ、不同 Δt。这些信息是 per-block 的(`getMaxLevel()`、per-block Omega)。单 block 只能有一种分辨率,**GPU 路径与 AMR 目前完全互斥**,而 AMR 是 FreeLB 的主打能力(README 首屏两张图)。 |
| **多 GPU** | 一块卡一个进程/设备时,block 是天然的分配与负载均衡单位。 |

### 2.2 不成立 —— 常见理由经核算站不住(诚实标注)

| 常见理由 | 核算结果 |
|---|---|
| **"单 block 装不下,显存放不下"** | 不成立。每 cell 约 93 B(POP 19×4 + rho 4 + u 12 + flag 1);6 GB 显存可放约 6.5e7 cell ≈ **400³**。单 block 能吃满整卡显存,显存不是分块的理由。 |
| **"复杂几何有大量 void,分块能省掉无效线程"** | **方向对,但不归多 block**。GPU 现在对 block 内全部 N 个 cell 起线程再按 flag 分派;省掉 void 的正解是**活跃索引压缩**(已有 `GenericvectorManager<std::size_t> BulkTaskIds(Geo.getBlockNum(), FlagFM, AABBFlag)` 这种 per-block 索引表),一个 block 也能做,launch 数仍是 1。这是与多 block **正交**的优化,应单独立项,不应算在多 block 的收益里。 |
| **"多 block 能提速均匀网格"** | 反作用。均匀密网格单卡下,多 block 只增加 launch、halo、overlap 三项纯支出(见 §3)。 |

### 2.3 结论

> **多 block 的目标场景是 AMR 与多 GPU,不是给均匀网格提速。**
> 立项时如果拿"提速"当目标,会在 M1 阶段就发现自己做出来的东西比现在慢 30%–100%,然后进退两难。先把目标锚定在"能力"上。

---

## 3. 性能损耗分析

### 3.1 实测(本次新做)

探针:`benchmarks/multiblock_probe/gpu_multiblock_probe.cu`(空 kernel 测纯 launch;DRAM-bound copy 测"固定总工作量拆成 N 次 launch"的真实代价;两者都做 CUDA Graph 对照)

```
MSYS_NO_PATHCONV=1 wsl.exe -- bash -lc \
  '/usr/local/cuda/bin/nvcc -O3 -arch=native \
   /mnt/d/WorkBuddy/wsl/gpu_multiblock_probe.cu -o /tmp/probe2 && /tmp/probe2'
```

**[A] 空 kernel,µs/step(500 reps + 20 warmup)**

| N launches | 逐次 launch | CUDA Graph | ns/launch | Graph 加速 |
|---|---|---|---|---|
| 1 | 7.68 | 7.73 | 7676 | 0.99x |
| 2 | 14.74 | 7.18 | 7370 | 2.05x |
| 4 | 31.01 | 7.50 | 7753 | 4.13x |
| 8 | 61.39 | 9.64 | 7674 | 6.37x |
| 16 | 122.54 | 16.16 | 7659 | 7.58x |
| 32 | 239.30 | 30.21 | 7478 | 7.92x |
| 64 | 487.54 | 60.11 | 7618 | 8.11x |
| 128 | 848.23 | 108.13 | 6627 | 7.84x |

边际成本:**逐次 6.6 µs/launch;Graph 0.79 µs/launch(差 8.4x)**。

> ⚠️ 7.6 µs/launch 这个绝对值偏高,WSL2 的 GPU 直通会放大 launch 延迟;原生 Linux 预期 2–3 µs。**绝对值请以原生 Linux 复测为准,但"边际线性 + Graph 能压 8x"这个结构性结论不依赖绝对值。**

**[B] DRAM-bound 总工作量固定(4M float,32 MB/step,308 GB/s 满速)拆成 N 次 launch**

| N launches | µs/step | GB/s | vs N=1 | Graph µs | Graph GB/s |
|---|---|---|---|---|---|
| 1 | 108.85 | 308.3 | 1.00x | 108.26 | 310.0 |
| 2 | 111.97 | 299.7 | 1.03x | 110.00 | 305.0 |
| 4 | 121.99 | 275.1 | 1.12x | 114.94 | 291.9 |
| 8 | 131.98 | 254.2 | 1.21x | 113.35 | 296.0 |
| 16 | 321.03 | 104.5 | **2.95x** | 121.99 | 276.9 |
| 32 | 575.43 | 58.3 | **5.29x** | 140.31 | 239.1 |
| 64 | 1334.98 | 25.1 | 12.26x | 175.78 | 190.9 |
| 128 | 986.70 | 34.0 | 9.06x | 245.52 | 136.7 |

读法:

- N ≤ 8 时损耗温和(1.2x 以内),**N ≥ 16 断崖式恶化**;
- 恶化幅度**大于** [A] 的纯 launch 叠加(16 launch: [A] 预测 +122 µs,实测 +212 µs),说明除 CPU 侧入队外还有 GPU 侧 per-kernel 开销(kernel 收尾/前端),**机制待 nsys 确认**,不要拿 [A] 直接外推;
- N ≥ 32 两次复跑离散度约 ±20%,该区间数字只作量级参考;
- **Graph 把 N=16 从 2.95x 拉回 1.12x、N=32 从 5.29x 拉回 1.29x** —— 这是本项目里唯一能把多 block 从"亏"变成"平"的手段。

### 3.2 损耗项清单

| 损耗项 | 量级 | 随什么变化 | 能否消除 |
|---|---|---|---|
| **kernel launch** | 实测边际 6.6 µs/launch(WSL);每 block 每步 2 次 | ∝ block 数 | **CUDA Graph(→0.79 µs)**;再把 stream 的 N 次 launch 合并成 1 次(grid = N×q),2N → N+1 |
| 波前量化(尾部) | `ceil(n/128)` CTA 向上取整到常驻 150 CTA 的整数倍;125k cell → 约 7% | 小块显著 | 靠 block 尺寸下限约束 |
| 占用率不足 | 常驻容量 = 30 SM × 5 CTA(96 regs × 128 thr = 12288 regs/CTA,65536/12288 = 5)= **150 CTA = 19200 cell**;block 小于此值时 SM 吃不饱 | 小块 | 同上 |
| **overlap 冗余** | `(B+2)³/B³ − 1`:50³ → 12.5%,25³ → 26%,16³ → 42% | 块越小越差 | 块尺寸下限;overlap 已是参数 |
| **halo 交换** | 流量 O(N²ᐟ³) vs 计算 O(N),量级上小;**但若走 D2H→CPU comm→H2D,每步 ~9 MB + 两次同步 ≈ ms 级,直接毁掉全部收益** | 面体比 | **必须写 device 端 halo kernel**,宿主 `normalCommunicate` 不能直接用 |
| RegPop 寄存器驻留 | 96 regs / LOCAL:0(D3Q19 FP32,已实测);运行时按 blockIdx 动态分派会打断驻留 | 实现相关 | 每 block 一个模板实例化的 kernel,block 索引作为**编译期**实参,不做运行时分派 |
| 显存/对象开销 | 每 block 约 5–6 次 `cudaMalloc` + H2D | ∝ block 数 | 一次性分配大块 + 指针偏移,初始化期成本,可接受 |

### 3.3 换算到真实 solver(100³,FP32 2211 MLUPS ⇒ 452 µs/步)

以 [B] 的**实测附加开销**线性叠加(开销与工作量无关),每 block 2 次 launch:

| block 数 | launch 数 | 附加 µs | 预计 µs/步 | 相对现状 | Graph 后 |
|---|---|---|---|---|---|
| 1 | 2 | +3 | 455 | 1.01x | 1.00x |
| 2 | 4 | +13 | 465 | 1.03x | 1.01x |
| 4 | 8 | +23 ~ +48 | 475–500 | 1.05–1.11x | 1.02x |
| 8 | 16 | +191 | 643 | **1.42x** | **1.03x** |
| 16 | 32 | +466 ~ +516 | ~930 | **2.06x** | **1.07x** |
| 32 | 64 | +647 ~ +1226 | 1100–1680 | 2.4–3.7x | 1.15x |

**推论**:

- 裸实现的多 block 在 **block 数 ≥ 8** 时就不划算;**Graph 化后 block 数 ≤ 32 仍可控制在 15% 以内**。
- 网格放大后这个结论会反转:400³ 时每步 ~29 ms,64 次 launch 的 0.4 ms 只占 1.4%。**所以损耗的本质是"每 block 的格子数",不是"block 数"** —— 判据应该是 `cells_per_block ≥ 5e4 ~ 1e5`,而不是限制 block 数。

### 3.4 由此得到的设计约束

> **每 block 建议下限 ≈ 5×10⁴ cell(约 37³)**。
> 下限由三条同时给出:占用率地板 19200 cell、overlap 冗余在 50³ 以下开始超过 12%、launch 开销占比要求单块工作量 ≥ 数百 µs。

---

## 4. 实施计划

### M0 安全护栏(约 0.5 天,零性能影响,建议无条件先做)

| 项 | 内容 | 文件 |
|---|---|---|
| M0-1 | `BlockLatticeManager::CuDev*` 在 `BlockLats.size() > 1` 时**直接报错退出**,不再静默只算 block 0 | `block_lattice.hh` |
| M0-2 | 两个 cell-dynamics kernel 加 `if (idx < N)` 边界检查(补上已被注释掉的 N 重载的调用) | `block_lattice.h` / `.hh` |
| M0-3 | checksum 覆盖所有 block(`getBlockLat(0)` → 遍历) | `cavity3d.cu`, `field_checksum.h` |

回退:全是加断言,无行为变更(单 block 路径 checksum 应保持逐位一致)。

### M1 单卡同分辨率多 block(约 2–3 天)

| 步骤 | 内容 | 文件 | 预计量 |
|---|---|---|---|
| S1 | 恢复被注释的多 block 循环(`CuDevStream` / `CuDevApplyCellDynamics` / `...Reg` 四处) | `block_lattice.hh` | ~15 行 |
| S2 | stream launch 合并:一个 kernel,grid = Nblock × q,替代每 block 一次 `<<<q,1>>>` | `block_lattice.h` / `.hh` | ~20 行 |
| S3 | device 端 halo kernel:把 `CommDirection` / send-recv 索引表 H2D,写 gather/scatter kernel;**先只做同分辨率(直接拷贝,不插值)** | 新增 `cuda_comm.h` | ~150 行 |
| S4 | 每 block 一次大块显存分配,替换 5–6 次小 `cudaMalloc` | `block_lattice.h` | ~30 行 |
| S5 | benchmark:`block_num` 参数化 + MLUPS 扫描 | `cavity3d_cu` | ~30 行 |

**验证**:见 §5。回退:`block_num = 1` 时与现状逐位一致。

### M2 损耗消解(约 1–2 天,依赖 M1 实测结果)

| 步骤 | 内容 |
|---|---|
| S6 | CUDA Graph 化:每步一次 `cudaGraphLaunch`(block 数不变时可复用;block 数变化需重建) |
| S7 | nsys 定位 [B] 中断崖式恶化的真实来源(CPU 入队 vs GPU 前端),据此决定要不要继续做 S8 |
| S8 | (可选)把 collision 与 stream 合并进同一 kernel,launch 数 N+1 → 1 |

**门槛**:M2 做到 "block 数 ≤ 32 时损耗 ≤ 15%" 即停手,不要追求更多。

### M3 跨分辨率 / AMR(二期,单独立项,本计划不做实施)

三件难事,建议先只做可行性判断:

1. **界面插值**:不同 Δx 之间 halo 需要时空插值(`popIntp` 目前只有宿主版);
2. **局部时间步**:细层一个 Δt 走 2 步、粗层走 1 步 —— 与现有"每步全场两个 kernel"的主循环结构冲突,主循环要重写;
3. **旋转 streaming 失效**:`CyclicArray::rotate` 依赖块内均匀平移,粗细界面处必须退化成真正的 gather/scatter。

### M4 多 GPU(远期)

一个 rank 一张卡,block 到设备的映射 + 现有 MPI 通信复用。依赖 M3 的结论。

---

## 5. 验证方案

**正确性(三级,与 FP16_PLAN 同一套判据)**

1. **回归**:`block_num = 1` 时,GPU 输出 checksum 与 M0 之前**逐位一致**;
2. **跨实现**:GPU 多 block vs CPU 多 block(同 `block_num`),1000 步、100³、float——目标 sum 相对差 ~1e-7 量级(单 block GPU vs CPU 已知 2.3e-7);
3. **物理不变性(最强判据)**:把同一个算例分别按 1 / 2 / 4 / 8 block 跑,收敛后比较速度剖面与总质量 —— **分块不应改变物理结果**,任何随 block 数漂移的量都说明 halo 有问题。

**性能**

4. `block_num ∈ {1,2,4,8,16,32,64}` × 网格 `{100³, 200³}` 的 MLUPS 曲线,给出"损耗 vs 每 block 格子数"曲线而非 vs block 数;
5. nsys 抓 launch gap、SM 占用率、halo kernel 占比;
6. Graph 前后对照(§3.1 [B] 已有对照方法)。

**稳定性**

7. 8 block 长跑 3000 步不 NaN、质量守恒;对照 FP16 的长跑失效边界记录方式。

---

## 6. 风险与回退

| 风险 | 说明 | 处置 |
|---|---|---|
| **静默错误**(最高) | `block_num>1` 只算一块、越界写、launch 失败未检查 —— 已发生过一次同类事故 | M0 无条件先做 |
| 通信走宿主 | 一旦偷懒用 D2H/H2D 做 halo,每步 ms 级,收益全无且不易察觉 | S3 是硬门槛,不做完不许进性能评测 |
| CUDA Graph 与动态性冲突 | block 数/网格变化需重建 graph;AMR 动态加密会频繁触发 | M2 限定"静态块结构"场景;动态加密留给 M3 评估 |
| **过度工程** | 若 M1 实测 block 数 ≥ 8 就亏且 Graph 也救不回来,应**止步于 M0 + 活跃索引压缩**,把多 block 只保留给 AMR 这一"必须有"的场景,不作为默认路径 | §7 门槛 |

---

## 7. 决策门槛(Go / No-Go)

**Go** —— 满足任一:

- 需要 **AMR on GPU**(FreeLB 的核心能力,目前 GPU 路径完全不支持);
- 需要 **多 GPU**;
- 复杂几何算例中,void 比例高到值得做几何贴合的块分解(**注意:先用活跃索引压缩试,那条路 launch 数不变**)。

**No-Go**:

- 只为"给均匀网格提速" —— 不做。FP16 已经拿到实测 1.86x(2211 → 4113 MLUPS),比多 block 靠谱一个量级,先做那个。

**M2 通过线**:block 数 ≤ 32、每 block ≥ 5×10⁴ cell 时,多 block 相对单 block 损耗 ≤ 15%。达不到就把多 block 定位为"AMR 专用路径",在文档里写明它对均匀网格是负收益。

---

## 附:探针

`benchmarks/multiblock_probe/gpu_multiblock_probe.cu` —— 空 kernel launch 开销、固定工作量拆分代价、CUDA Graph 对照,三者合一。§3.1 全部数据由它产出,可直接复跑。

---

## 8. 补充:复杂几何下的 void 代价与显存取舍(2026-09-26)

针对"单卡、短期内不上多 GPU、但复杂几何包围盒远大于实体"的具体纠结。结论是
**先别在多 block 上做取舍 —— 有一笔账还没算,而且它不是显存账。**

### 8.1 void 单元当前是全额付费的,付的是带宽

默认路径(`--reg`,即 `cudev::Cell<..., cudev::RegPop>`):

```cpp
cudev::Cell<..., cudev::RegPop> cell(idx, blocklat);  // 构造函数里:
                                                      //   getPopArray 解析一次地址
                                                      //   PopCache::load: v[d] = p[d][0]
                                                      //   → 19 次 load 已经发生
CELLDYNAMICS::Execute(flagarr->operator[](idx), cell);   // flag 分派在这之后
cell.flush();                                         // → PopCache::store,19 次 store 也必然发生
```

`src/data_struct/cell.h` 里 RegPop 特化下的构造函数无条件搬运全部 q 个分布函数
(这是寄存器驻留的设计前提,不是 bug);`flag` 判定在搬运之后。

对照 baseline 路径(`--base`):`cudev::Cell::operator[]` 只返回全局内存引用,
惰性求值,void cell 不触发任何访存。

**推论**:

- `--base`:void 单元不产生流量;
- `--reg`(默认):void 单元产生完整的 152 B/cell 流量;
- 100³ cavity 是 100% 填满的,所以这笔开销在现有基准里**完全不可见**;
- 按有用格点算,**有效 MLUPS = 实测 MLUPS × φ**(φ = 实体占包围盒的体积比)。
  φ = 0.3 时,2211 MLUPS 实际只有约 663 有效 MLUPS —— 已经是 3x 量级的隐形损失,
  比多 block 的 −30% 大一个数量级。

### 8.2 于是正确的顺序是:先拿纯收益的那一步

**活跃索引压缩 + 延迟搬运**(约 1 天,与多 block 正交;即让 RegPop 的 load 延后到 flag 判定之后)

- kernel 不再从 `0..N−1` 起线程,而是从活跃 cell 索引表起线程
  (已有 `GenericvectorManager<std::size_t> BulkTaskIds(Geo.getBlockNum(), FlagFM, AABBFlag)`,
  本来就是 per-block 的);
- **旋转 streaming 让这件事异常干净**:pull 模式下每个线程只访问自己那条 cell 的
  19 个 slot(`data_d[i]`),不 gather 邻居。所以压缩索引**不需要邻居表、不需要改
  CyclicArray、不需要 halo**;
- 收益:每步流量 ∝ N_active;代价:仍是 1 次 launch,访存局部性取决于活跃格点在
  最快轴上的连续游程长度;
- **不省显存**。

这一步在"做多 block"和"不做多 block"两种结局下都赢,不可能亏。

### 8.3 显存这笔账 —— 以及明确的门槛

每 cell 约 93 B(POP 19×4 + rho 4 + u 12 + flag 1),6 GB ⇒ 约 6.5e7 cell ≈ 400³。

```
N_bbox × 93 B ≤ 0.7 × 可用显存   →   只做 §8.2,不要碰多 block(维持 §7 的 No-Go)
否则                              →   显存成为真瓶颈,才谈 tiled 分配
```

理由很直白:**浪费你已经拥有的显存没有成本。** 只有当"算不了"或"必须降分辨率"
时,tile 才值得 —— 降分辨率损失的精度远大于 halo 的 20–40%。

### 8.4 若显存真卡住:tile 的尺寸张力是真的

tile 边长 m、实体体积 V、表面积 A(格点数):

- 需分配显存 ≈ V·(1 + 6/m) + A·m
- 每步时间 penalty ≈ 6/m(halo 冗余计算)+ launch 开销
- VRAM 最优解 m* = √(6V/A):球(R=100)→ m ≈ 14;厚 8 格的薄板 → m ≈ 5

**VRAM 最优点落在小 tile(halo penalty 40%–100%),而时间 penalty 要求大 tile。**
这是真实且无法消除的张力 —— 用户原本的顾虑并非错觉。

**但张力里的 launch 项不是本质的,是实现决定的。**破局:不要写成"每 block 一次
launch",而写成**一个 kernel 遍历 tile 表**(grid-stride over tiles),launch 数
恒为 O(1)/步,§3.1 实测的 6.6 µs/launch 就不再是 tile 数的函数。剩下的 6/m
躲不掉,所以 m 的选取准则变为:**在显存预算允许范围内取最大的 m**,而非 VRAM 最优的 m*。

### 8.5 第三条路已被实测否掉:虚存稀疏提交

设想:保持稠密索引不变(所有常偏移、rotate streaming、合并访存一概不变),只给非空
VA 区间提交物理页 —— 理论上是零性能代价的省显存。仓库里已有现成 VMM 代码:
`StreamMapArray::setDevMap()` 的 `cuMemAddressReserve` + `cuMemCreate` +
`cuMemMap` + `cuMemSetAccess`。

**实测否定**(本卡,`benchmarks/multiblock_probe/vmm_granularity_probe.cu`):

```
device=NVIDIA GeForce RTX 3060 Laptop GPU  VMM supported=1
granularity MIN=2097152 (2.00 MB)  REC=2097152 (2.00 MB)
FP32 POP: 524288 cells/chunk → 等效立方边长 80.6;nx=400 时是 3.27 个 xy 平面
FP16 POP: 1048576 cells/chunk
```

单个 chunk 覆盖单个方向数组上 524288 个**连续** cell。在线性 (i,j,k) 序下,只要几何
沿最慢轴铺开,一个 chunk 都省不掉。要拿到细粒度必须改成 Morton 序,而那会破坏常偏移
→ 旋转 streaming 失效 → 退回 gather streaming。**在不改 streaming 机制的前提下这条路
不成立。**

副作用知识:`getPageAlignedCount` 会把 19 个方向各自向上取整到 2 MB,这是
StreamMapArray 在小算例上会多占一大块显存的原因 —— 也是当初默认切到 CyclicArray
的一个合理理由。

### 8.6 顺带:最便宜的显存杠杆已经在手

FP16 存储已落地(§见 docs/FP16_PLAN.md,4113 MLUPS):每 cell 93 B → 约 55 B,
**可算规模直接放大 1.69 倍**,且不引入任何 tile/halo/launch 复杂度。
在动多 block 之前,先把这个杠杆用满。

### 8.7 决定这件事的两个测量(都很快)

1. **判定 8.1 是否成立**:现有代码不动,cavity3d 里把一部分 cell 的 flag 设为 void
   (保持 N 不变),看 MLUPS 是否变化。**若 MLUPS 基本不动 → 证实 void 单元在全额
   消耗带宽**,那么 §8.2 就是当下收益最大的一件事。
2. **拿到真实几何的 φ 与 A/V**:对实际复杂几何统计 `getTotalCellNum()` vs 按 flag
   计数的 active 数。φ > 0.5 且 400³ 以内装得下 ⇒ §8.2 做完即止,多 block 永远不做。

在这两个数字出来之前,关于多 block 的所有性能讨论都只是估值。
