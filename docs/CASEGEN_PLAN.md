# 算例生成器与算子 DSL（tools/casegen）的决策记录

状态：Phase 1–4 已落地（Phase 4 为 geometry file support）；本文件只记录决策、被否决方案与延迟范围，用法与字段表见 `tools/casegen/README.md`｜日期：2026-11-06｜目标：回答“为什么是现在这个形状”，不重复用户手册

> 退出条件：`tools/casegen/README.md` 与代码能独立解释用法和实现时，本文件只剩“为什么”。
> 改实现不需要改本文件；改变下面某个决策才需要。

## 1. 边界

生成器只负责装配：字段集合、几何涂色、边界、算子选择、主循环骨架、`.ini` 与 `Makefile`。
物理内核来自内置目录或 `.op` 算子。产物是一次性脚手架，用户拥有并可以继续手改。
除 §2 记录的 periodic comm 存储修复外，不改 `src/`；不改 `make.mk`、`tests/run_tests.sh`；
不做 CA、自由表面、AMR。

这是有意的窄口：结构是编译期类型装配，运行期参数已经在 `.ini`；把两者都搬进一门新语言，
收益是省几百行誊写，代价是重新实现类型检查、诊断和一层永久维护面。

## 2. 关键决策与证据

- **Python 直接发射 C++，不复用 csegen。** `csegen` 的输入 basename 硬编码为
  `moment.h/equilibrium.h/force.h`（`plugins/freelb/ur_emit.cpp:116-127`），先泛化它会把
  submodule 的 CI 与行号契约一起拖进来。算子体是逐格算术，C++ 编译器已经能做常量折叠与
  循环展开；DSL 先保证语义与 `file:line` 诊断，CSE 留作后续优化。
- **两套语法，按名字链接，不合成一个文件。** 结构 deck 是类型/列表/顺序，内核是表达式；
  `third_party/cse` 前端只有 7 种语句，顶层未建模 token 会被原样放过
  （`src/frontend/parser.h:33-42`），一个 `[section]` 就能让整段静默跳过。`.case` 只写
  `[operators] name = path`，`.op` 独立解析、独立报错。
- **`box` 是 AABB 端点，不是格号闭区间。** 格心在 `(i+0.5)*cell_len`，
  `AABB::isInside` 是闭区间 + eps（`src/geometry/basic_geometry.hh:41-49`）。`0,0,32,32`
  覆盖格号 `0..31`，与 `cavity2d` 的 `AABB{0, Ni*Cell_Len}` 一致；按格号闭区间实现会多涂
  一排，1000 步基线必然失败。
- **`boundary` 只允许全域。** `SetupBoundary` 的判据是“任一邻居在 box 外就涂”，且不检查
  格子自身（`src/geometry/block_geometry2d.hh:55-76`）；非全域会静默涂到 box 外。
- **一次性脚手架，不做“配置即源、每次重生成”。** 用户手改一次就分叉，还会把生成代码与
  `src/` API 的漂移带进构建。产物归用户；`--check` 只用于 golden 目录防漂移。
- **校验和要锁循环相位，不只锁聚合量。** POP 的 sum/sumsq 对置换不变，单靠它不能区分
  `Stream→collide` 与 `collide→Stream`；Phase 1 因此同时钉 RHO sumsq。CUDA 沿用
  `cavity2d_cu` 的 `collide→Stream`，与 CPU 的 `Stream→collide` 相比边界相位差使
  `maxabs` 差约 `4.9e-4`，但 RHO sumsq 只差 `1.5e-7`，所以跨后端只钉 RHO。
- **周期边界的 `BasicComm&` 必须稳定。** `LatticeComm` 持有 `BasicComm&`；
  `BasicCommSet` 用 `std::deque` 存储，`emplace_back` 不会使已有元素引用失效，
  这是 casegen 侧唯一的 `src/` 改动；`reserve(16)` 那类容量补丁不再需要。
- **Geometry 是编译产物，不是 case 目录里的 opaque binary。** `.case` 的
  `[geometry.source.*]` 是唯一源描述；casegen 把它渲染成可 review 的 `.geo` 文本，
  `tools/geometry_build` 编译成 `.flbgeo`，生成的 Makefile 在链接前重建。`--check` 比较
  `.geo`，`*.flbgeo` 被 `.gitignore` 忽略；prebuilt `file` 模式保留，用于外部工具生成
  或暂不纳入 recipe 的 mask。
- **Source 顺序是格式语义。** `.geo` 的 `[source.<name>]` 按文件顺序应用，后写的 solid
  可以覆盖先写的 domain；解析器因此不能排序，也不能换成不保序的 `.ini` reader。
- **Geometry 的 patch/entity 尚未接 task/BC（B4）。** `.flbgeo` 持久化
  role/entity/material、boundary patch 与 surface band，但生成的 case 仍从 `[tasks]` 和
  `[region.<name>]` 选任务；patch 目前只是数据，不参与边界条件。

## 3. 被否决的方案

- **一个文件同时写结构与内核**：两套语法、行号与静默跳过，见 §2。
- **运行期解释器**：逐格分派要求内联、`__any__`、无分配，CUDA kernel 要求同源。
- **配置即源、每次重新生成**：用户手改即分叉；见 §2。
- **复用 `csegen` 作为算子后端**：basename 硬编码，submodule 契约；见 §2。
- **固定 driver + 生成头**：cavity2d/3d 的 diff 只有约 40 行，为它引入 hook 体系不划算。
- **提交 `.flbgeo` fixture 作为真源**：二进制、native-endian、scalar-tagged，review 不可读，
  格式或标量类型升级要重建所有 fixture。改为提交 `.geo` 文本，由工具重建。
- **casegen 直接写 `.flbgeo`**：会把 rasterise、band、patch 的语义复制到 Python 侧，与 C++
  builder 形成两份实现。`.geo` + `tools/geometry_build` 保持单一实现。

## 4. 延迟与未实现

- `open` / `anti_pressure`：生成期明确报错，`--list` 不列出。退出条件：有真实开边界用例
  与固定基线后再接入。
- `plain` + MPI：生成代码在 `mpi().getSize() > 1` 时明确报错，直到有 plain 的 `mpirun`
  校验和测试。`periodic` 的多层/多轴：未验证。
- FP16：未验证。按 `AGENTS.md` §3 的 CUDA/FP16 行手工跑过后再声明。
- **Geometry 消费端**：patch/entity 接 task/BC（B4）与最小壁厚校验未做；2D 路径保持
  原样。`.vti`、CUDA device mirror、`examples/` 迁移仍未做。

## 5. 验收基线

- `python3 -m unittest discover -s tools/casegen/tests`：185 项。
- 13 个 golden 目录：
  `tests/feature/casegen{,_op,_dsl,_dsl_math,_3d,_3d_cu,_cu,_force,_mrt,_les,_mpi,_plain,_geo_box3d}`，
  每个 `make` 前跑 `--check --selftest`。
- `make test-casegen`：自动发现全部 golden；CUDA 按 nvcc/GPU 跳过，MPI 目录跑
  `mpirun -np 2`。本机 `PASS 14 / FAIL 0 / SKIP 0`（`check.sh` 把单测汇总也计为一项）。
- `make test-cpu`：PASS 36 / FAIL 0 / SKIP 4；`tests/regression/cavity2d` 仍为
  `7728.702031`。
- 固定值表在 `tools/casegen/README.md` 的 “Exact checksums”；实现细节以代码为准。

## 6. 维护

- 改实现（加一个算子、改一段发射）：只改代码、README 与测试，不改本文件。
- 改决策（例如真的引入运行期解释器、把 DSL 合进 `.case`）：改本文件，并写明新的证据。
- 新目录算子必须有 golden 或一次手工验证记录，否则不进 `--list`。
- 改 `.geo` 字段或 `.flbgeo` 的版本/标量标签属于格式变更：先定兼容策略与 fixture 重建
  路径，再改代码。

## 7. 后续方向

优先级按“收益 / 维护成本”排；任何新功能必须有 golden 或一次手工验证记录，否则不进 `--list`。

| 方向 | 价值 | 代价 | 判断 |
| --- | --- | --- | --- |
| CI 接入 + `make test-casegen` | 高：golden 是唯一防漂移手段 | 0.5–1 天 | 已落地 |
| DSL 数学函数 / 向量分量 / 白名单字段 | 高：解锁真实自定义算子 | 1–1.5 周 | 已落地 |
| 3D CUDA 等价基线、plain+MPI 取舍 | 中：补两个验证洞 | 1–2 天 | 已落地 |
| 生成 Makefile 的 `-O3 -march=native` 选项 | 低中：避免生成算例比手写慢 | 0.5 天 | 已落地 |
| `open` / `anti_pressure` | 中高：最大的缺失边界族 | 1–2 周 | 有真实开边界用例再做 |
| geometry patch/entity → task/BC（B4） | 高：让 patch 真正参与边界 | 依赖 link-BC 与 B4 基线 | 先于 `.vti` |
| `.vti` 导出 | 中：可视化/调试 geometry mask | 1–2 天 | 有消费者再做 |
| 稀疏几何的填充式分解（B3） | 高：避免整域块浪费 | 预算见 `COMPLEX_GEOMETRY_PLAN` §7 | 先于 B5 |
| 区域初值 `set_value` | 中：AABB 之外最实用的缺口 | 3–5 天 | 按需求 |
| 目录扩展（力变体 / 非牛顿 / 热） | 中：每个 1–2 天 + 一个 golden | 需求驱动 | 按用例加 |
| CA / 自由面 / AMR / 运行期解释器 | — | 数月 | 不做 |

### 已落地

- **自动化**：`check.sh` 自动发现 `tests/feature/casegen*`，断言单测数大于 0 且每个
  Makefile 保留 `--check --selftest`；CUDA 按 nvcc/GPU 跳过，MPI 目录跑 `mpirun -np 2`。
  CI 有 `casegen` 与 `cpu-tests` 两个 job。
- **DSL**：`sqrt/abs/min/max/clamp`、`u[k]`、`field(RHO/VELOCITY/FORCE/SMAGORINSKY)`；
  `force()` 是 `field(FORCE)` 别名。下标按 `q/dim` 检查。验收算例 `casegen_dsl_math`
  （POP `586.5034095`、RHO `1155.987548`）。
- **验证洞**：3D CUDA 基线 CPU(float) RHO `5831.932013`，GPU `5831.9294`（`4.5e-7`）；
  MPI 内点校验和 `casegen_mpi`（np1 = np2：POP `445.51566903028549`、RHO
  `1023.9878563745673`）；`plain` 多 rank 明确报错；`thread_num != 1` 生成期拒绝。
- **审计修复**：自定义 `.op` 进入 rho/U pass（`casegen_op` 与目录 BGK 逐位一致）；
  `f[9]`/2D `u[2]` 生成期报错；`cuda = 1` + `SMAGORINSKY` 生成期拒绝。
- **构建**：普通生成的 CPU Makefile 用 `-O3 -march=native`；golden 用 `-O2` 保持基线
  可复现。
- **Geometry file（Phase 4）**：`[geometry]` 支持 `file`、`domain_flag`、`void_flag`、
  `band_width`、`origin`；`[geometry.source.*]` 支持 box/sphere/stl，按文件顺序应用。
  casegen 生成 `.geo`，生成的 Makefile 调 `tools/geometry_build` 重建 `.flbgeo`。golden
  `casegen_geo_box3d` 为 16³、100 步；prebuilt `file` 模式保留。

### 暂缓

- `open`/`anti_pressure`：需要进/出口区域、压力边界和另一套主循环；没有真实用例时
  只增加组合矩阵。
- STL/CSG 直接导入：casegen 仍只生成 `.geo`，由 `tools/geometry_build` 栅格化；直接
  导入需要几何内核/布尔运算，暂不做。
- 目录扩展：每个条目都要一个 golden，没有对应算例时维护成本大于收益。
- CA / 自由面 / AMR：另一套骨架，不是 casegen 的增量。
