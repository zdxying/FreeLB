# FreeLB 的 CSE 代码生成器（`.ur.h`）

> **迁移说明**：本目录旧的解释器/优化器（约 2900 行）已被独立仓库 `cse`
> （`third_party/cse` submodule，`main` 分支）中的通用 DAG-CSE 引擎取代。
> 原 `freelb-port` 工作已合入 `main`。
> 本目录现在只做两件事：把引擎接进 FreeLB 构建、提供 `csegen <in.h> <out.h>`
> 驱动与数值验证脚本。引擎侧的完整状态见
> `third_party/cse/docs/freelb_port_status.md`。

## 概述

`csegen` 读取 `src/lbm/{moment,equilibrium,force}.h` 中带 `// @cse` 标记的模板
结构体（其 `apply(...)` 方法体内含按格子集展开的循环），对每个格子集生成完全
展开、公共子表达式已消除的 `.ur.h` 特化。生成结果与手写 `.ur.h` 在数值上一致
（`make verify`），在 `-D_UNROLLFOR` 下经影子包含目录替换手写版本。

支持 6 个格子集：D2Q5、D2Q9、D3Q7、D3Q15、D3Q19、D3Q27（`cs2 = 1/3`）。

## 架构总览

```
FreeLB src/lbm/*.h  (// @cse 标记)
        │
        ▼
third_party/cse 引擎
  ├─ region_extractor   提取 @cse 区域
  ├─ lexer / parser     → AST（函数、结构体、using、namespace）
  ├─ ur_emit            按结构体模板形参分类，逐 latset/逐 d 实例化
  │     每个实例走一遍：
  │     ├─ IRBuilder        AST → DAG IR（含向量降级、latset 常量折叠）
  │     ├─ PassManager      loop_unroll → lattice_resolve → constant_fold →
  │     │                   algebraic_simplify → [counter_prop] →
  │     │                   [reassociate] → cse → [recombine] →
  │     │                   value_prop → dce
  │     └─ CodeGen          发射方法体
  └─ .ur.h 组装（文件头 / #ifdef _UNROLLFOR / namespace / 偏特化）
        │
        ▼
generated/lbm/*.ur.h  (由 make.mk 影子包含)
```

`src/` 是通用引擎，提供静态/动态库 `libcse.a/.so`；`bin/cse`（对 `//@cse`
区域做 CSE，写 `<in>.cse`，见 `tests/`）由插件驱动
`plugins/freelb/cse_main.cpp` 构建，默认采用 FreeLB 配置。

## `.ur.h` 生成（`plugins/freelb/ur_emit.*`）

1. **标记与区域**：`region_extractor` 识别 `// @cse`（兼容 `//@cse`）；仅处理
   紧跟其后的一个结构体定义。
2. **结构体分类**（`classifyStruct`）：
   | 模板形参形态 | 类别 | 实例化方式 |
   |---|---|---|
   | `<typename CELLTYPE, ...>` | `Cell` | 单次（保留模板形参） |
   | `<typename T, typename LatSet>` | `TLatSet` | 每个 latset 一次 |
   | `<typename T, typename LatSet, 非类型 d>` | `TLatSetD` | 每个 latset × `d=0..dim-1` |
   | `<typename CELL, ...>` | `CellType` | 使用方（如 moment 的 11 个结构体） |
   | 其它 | `Unsupported` | 跳过并告警 |
3. **每实例管线**：构造 `IRModule` → `IRBuilder` → `PassManager::createDefault`
   （插件传入 `createLatticeResolvePass("LatSet", <set>)`，并在 `lowerVectors`
   时把 `createCounterPropPass(vectorLocalName)` 作为 `postAlgebraPass`）→
   `CodeGen::generateBody`。
4. **组装输出**：写 `#pragma once`、include、`#ifdef _UNROLLFOR`、`namespace`、
   `CELL<T, LatSet<T>, TypePack>` 偏特化，以及每个 latset（`TLatSetD` 再乘每个
   `d`）的方法体。

### 向量与 latset 常量

引擎核心（`src/`）不含 FreeLB 语义；以下行为均由 `plugins/freelb/config.h`
经 `CSEConfig` 钩子注入：

- `CSEConfig.lowerVectors` + `vectorDim`/`isVectorType`/`isVectorProducingCall`：
  把 `Vector<T, LatSet::d>` 局部变量按分量降级；分量式 `+ - * /`、点积、常量
  分量索引、方向向量 `latset::c<LatSet>(k)` 在 `lattice_resolve` 后变为常量。
- `CSEConfig.resolveName`（`freelb::resolveLatsetConst`）把
  `<LatSet>::q/d/cs2/InvCs2/InvCs4` 折为常量；`lattice_resolve` 把
  `latset::c/w/opp` 解析为常量，权重规范化为 `latset::w<LatSet>(k)` 符号，
  不在折叠中丢失。
- `counter_prop`（仅在 `lowerVectors` 时由插件作为 `postAlgebraPass` 启用）：
  解析直线计数器（`tensor[i] → tensor[0]…`）、降级向量局部索引
  （`unew[1] → unew_1`，经 `CSEConfig.vectorLocalName`）、折叠常量 `if`。
- 非类型模板参数经 `CSEConfig.constBindings` 绑定为常量。

## 引擎内部要点

- **前端**：`src/frontend/{lexer,parser,ast}`；`parser` 支持 `T{}` 值初始化、
  `x.template f<...>()`、`if constexpr`、模板/限定类型局部声明。`CSEConfig`
  （`src/frontend/cse_config.h`）控制假设（交换/结合/FP 重结合/别名）并提供
  通用项目钩子（`resolveName`、向量降级、`vectorLocalName`）。
- **IR**：`src/ir/`。DAG 节点经 FNV-1a 强哈希去重（`ir_module`）；`IRBuilder`
  完成 AST→IR 与向量降级；`ir_utils.h` 提供 `foldConst`/`substitute`。
- **Passes**：`src/passes/`。`pass_manager` 的默认顺序见上；`loop_unroll`
  支持非零起点（三角循环）并逐迭代重命名可变/不纯局部；`cse_pass` 的候选
  选择按（语句数, 节点 id）排序，保证输出与堆布局无关。
- **后端**：`src/backend/codegen`。`if constexpr` 单语句不加花括号，贴合
  FreeLB 风格。
- **插件**：`plugins/freelb/`。`config.h`（激进配置 + 纯函数判定 + latset
  常量/向量触发器钩子）、`lattice_resolve`、`ur_emit`、`cse_main`；
  `cuda_skip` 用于跳过 CUDA-only 代码。
- **测试**：引擎 `make test` 运行代价回归、数值验证、`csegen` 冒烟，并在存在
  FreeLB checkout 时运行本目录的 `verify_*.py` 与
  `tests/verify/check_lattice.py`（latset 表防漂移）。

## FreeLB 集成

- `third_party/cse`：引擎 submodule（`main` 分支），远端
  `git@github.com:zdxying/cse.git`。
  更新到最新：`git submodule update --remote third_party/cse`。
- `make.mk`：当编译选项含 `-D_UNROLLFOR` 时，依据 `UR_CSE_BASES`
  （默认 `lbm/moment lbm/equilibrium lbm/force`）生成到 `generated/`，并把
  `generated/` 作为影子包含目录；未列出的头仍使用 `src/lbm/` 手写 `.ur.h`。
- 本目录目标：
  | 目标 | 作用 |
  |---|---|
  | `make` | 构建引擎并复制 `csegen` |
  | `make gen` | 由 `src/lbm/*.h` 生成 `out/*.ur.h` |
  | `make verify` | 与 `reference/*.ur.h` 数值对比（三个头 × 6 latset） |
  | `make gen-refs` | 用当前 `src/lbm/*.ur.h` 刷新参考快照 |
  | `make install` | `verify` 通过后把生成结果覆盖 `src/lbm/*.ur.h` |
- `reference/` 固定了手写 `.ur.h` 快照，使 `install` 之后 `verify` 仍具参考意义。
- `verify_*.py` 是数值验证契约（不可随意修改）：解析生成代码并对手写参考做
  数值比较。

## 文件结构

```
tools/cse/
  Makefile            引擎构建/接线、gen/verify/install
  csegen              由引擎复制而来的驱动二进制（构建产物）
  verify_equilibrium.py
  verify_force.py
  verify_moment.py
  reference/*.ur.h    手写参考快照
  out/*.ur.h          gen 输出
  PORT_STATUS.md      指针：迁移状态的唯一权威副本在引擎侧
  DESIGN.md           本文档
```

引擎（`third_party/cse`）——`src/` 无 FreeLB 语义，FreeLB 全部在
`plugins/freelb/`：

```
src/frontend/   lexer parser ast cse_config region_extractor
src/ir/         ir_module ir_builder statement ir_utils dag_node
src/passes/     loop_unroll counter_prop constant_fold
                algebraic_simplify reassociate cse
                expr_recomb value_prop dce pass_manager
src/backend/    codegen
src/analysis/   cost_model
plugins/freelb/ config lattice_resolve ur_emit ur_emit_main cse_main cuda_skip
tests/          run_tests.sh
tests/fixtures/ basic_cse features namespace_case
                equilibrium_d3q19 safety_cases
tests/verify/   verify_equilibrium verify_safety check_lattice.py
tests/csegen/   equilibrium force moment
```

## 使用方式

```bash
# 生成并验证（不修改 src/lbm）
cd tools/cse && make verify

# 生成到 generated/ 供 -D_UNROLLFOR 示例使用
cd tools/cse && make gen
cd ../../examples/cavity3d && make        # FLAGS 需含 -D_UNROLLFOR

# 用生成版覆盖手写 .ur.h
cd tools/cse && make install
```

## 相关文档

- `third_party/cse/docs/freelb_port_status.md`：迁移完成项与 TODO（唯一权威）。
- `third_party/cse/docs/architecture.md`：引擎通用架构。
- `PORT_STATUS.md`：本目录内的占位指针，指向上方引擎文档。
