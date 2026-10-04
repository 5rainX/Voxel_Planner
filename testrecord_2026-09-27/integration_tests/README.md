# 巴片规划集成测试与原始样本生成

这套脚本驱动师姐的原始 C++ 主程序，由它调用你提供的 Voxel_Planner，再执行师姐的细规划和现有约束检查。两份原项目源码均未修改。

调用链：

```text
run_tests.py
  -> hybrid_validation_cli.exe（师姐 src/main.cpp）
     -> voxel_planner::loadMap / findPaths（你的库）
     -> BusbarPlanner::plan（师姐的细规划和内部检查）
     -> TXT / CSV / OBJ
  -> results.jsonl + accepted.jsonl + summary.csv + report.md
```

## 1. 已验证的范围

- 5 个固定案例 + 10 个固定随机种子的合成案例，共 15 次真实 C++ 调用。
- 13 次通过细规划，得到 13 组粗路径—细路径配对；封闭墙与起点占据两个负例未进入候选数据集。
- 12 个 Python 回归测试通过，覆盖超时、崩溃、缺失输出、粗规划异常被 NO_PATH 覆盖、错误输入和禁止回退等情况。
- 原 Voxel_Planner 的 `test_sdf`、`test_elastic_band` 两项 CTest 通过。
- 本机验证环境：Python 3.9.2、MSVC 14.51、C++17、Ninja、OpenMP LLVM 模式。

详细实际结果见 `results/report.md`，数据见 `results/accepted.jsonl`。

**FINE_ACCEPTED 只表示当前师姐细规划器接受了结果。** 这不是全局最优证明，也没有增加独立 CAD、连续全尺寸碰撞或企业工艺认证。只有明确了你们对“标准解”的定义，才适合把这里的候选样本提升为正式标准数据。

## 2. 直接运行已经编译好的测试程序

在本目录打开 PowerShell，执行：

```powershell
python .\run_tests.py --exe .\bin\hybrid_validation_cli.exe --manifest .\cases\cases.json --output .\my_smoke_run --timeout 30 --fine-timeout 5 --threads 2 --export-obj
```

`--output` 必须是新目录，防止读到旧的成功文件。脚本不覆盖历史实验。

本机已安装该程序所需的 VC++ 和 OpenMP 运行库；换到其他机器可能需要相应运行库或重新编译。发生 DLL 缺失时会记为 `PROCESS_ERROR`，不能解释为规划无解。

## 3. 随机生成场景并积累样本

```powershell
python .\make_cases.py --output .\new_cases --random-cases 100 --seed 20260927
python .\run_tests.py --exe .\bin\hybrid_validation_cli.exe --manifest .\new_cases\cases.json --output .\new_results --timeout 60 --fine-timeout 10 --threads 2
```

这会生成 5 个固定案例和 100 个随机案例。地图为 48×32×24 的小型体素场景；宽厚为 2 和 1，间隙为 0，参数仅用于软件集成测试，不是企业制造参数。

随机生成器目前变化障碍布局与端点位置，宽厚和端部姿态保持固定；它不代表覆盖真实电池包分布，也不代表已经获得适合正式训练的数据规模。

`--export-obj` 可选。启用后，每个案例目录保存粗路径、粗路径截面形状和成功细路径的 OBJ，可用 MeshLab 等工具打开。

## 4. 接入真实地图

复制 `cases/cases.json`，保留需要的条目后修改：

```json
{
  "schema_version": 1,
  "units": "voxel",
  "cases": [
    {
      "id": "real_case_001",
      "map": "map.txt",
      "start": [8, 16, 12],
      "goal": [40, 16, 12],
      "args": [
        "--width", "2", "--thickness", "1", "--clearance", "0",
        "--min-begin-length", "1", "--min-end-length", "1",
        "--flat-bend-factor", "1", "--vertical-bend-factor", "1",
        "--twist-factor", "1", "--bend-factor", "1",
        "--start-normal", "0", "0", "1", "--end-normal", "0", "0", "1",
        "--start-tangent", "1", "0", "0", "--end-tangent", "1", "0", "0",
        "--no-bend45", "--heuristic", "3"
      ]
    }
  ]
}
```

这些坐标和参数只是格式示例，请替换成真实输入。`map` 相对于清单所在目录解析，也可使用绝对路径；`id` 只允许英文字母、数字、下划线和连字符。

地图格式必须是两份代码共同支持的稀疏 TXT，UTF-8 无 BOM：

```text
map: 48 32 24
24 10 5 0
24 11 5 0
```

第一行是网格尺寸，后续列出占据体素。师姐代码中 flag=0 为实体障碍、flag=1 为间隙碰撞体素；当前 Voxel_Planner 将每一条记录都视为占据声明，忽略末尾 flag 的细分语义。本测试保留原始输入，不擅自转换。

师姐主程序不能直接加载 VTK，因此集成测试使用 TXT。脚本会先检查格式与逐坐标边界，再复制为案例目录中的 `input.txt`；原始地图不会修改。

`args` 是师姐主程序的参数数组，支持其 `--via`、折弯 case 参数、候选数量等已有选项。脚本强制开启 `--voxel-planner-release`，未显式指定候选切向组时增加 `--endpoint-tangent-candidate 0 0`。这样粗路径失败时不会静默退回不使用你们库的旧流程。

会绕过或更换调用链的 `--coarse-only`、`--no-release-initial-path`、`--astar-jump` 等选项被拒绝。端点旋转请使用成对的 `--endpoint-tangent-candidate`，不要混用旧的单端旋转参数。

## 5. 如何看结果

| 状态 | 含义 | 是否进入 accepted.jsonl |
|---|---|---|
| FINE_ACCEPTED | 有实际粗路径，细规划返回成功，对应导出完整 | 是 |
| COARSE_NO_PATH | 本次粗规划未返回候选 | 否 |
| FINE_NO_PATH | 有粗路径，但细规划未返回成功结果 | 否 |
| FINE_TIMEOUT | 细规划报告 PLANNER_TIMEOUT | 否 |
| PROCESS_TIMEOUT | 整个原生进程超过外层时限 | 否 |
| COARSE_ERROR | 粗规划异常或预算耗尽等错误，保留原始原因 | 否 |
| PROCESS_ERROR | 原生进程启动/退出异常，例如 DLL 缺失 | 否 |
| INPUT_ERROR | 地图、坐标或请求参数格式有误 | 否 |
| OUTPUT_ERROR | 返回码、状态文件或导出内容不一致/不完整 | 否 |

`COARSE_NO_PATH`、`FINE_NO_PATH` 是算法在本次参数与预算下的搜索结果，不是数学不可行证明。

固定案例使用 `expected` 做回归检查。`NO_SOLUTION` 表示允许在粗或细阶段拒绝，但不把超时和异常算成正确拒绝。随机案例没有预设答案。

批处理退出码：0 表示已声明的预期均满足且没有运行/输出错误，**不表示所有案例都找到路径**。退出码 1 表示预期不符或出现超时、输入、运行、输出错误；具体看报告。

每次运行包含：

- `run.json`：二进制哈希、构建来源、线程和时间预算。
- `manifest.json`：本次输入清单。
- `results.jsonl`：所有尝试，包括失败；每行一个案例。
- `accepted.jsonl`：每行一个成功候选的配对，可能一场景有多行。
- `summary.csv` / `report.md`：汇总表。
- 每案例目录：输入副本、实际命令数组、stdout/stderr、候选表、路径/姿态/分段信息、可选 OBJ。

成功样本字段包括 `coarse_centerline`、`coarse_poses`（nx,ny,nz,tx,ty,tz）、`fine_segments`、原请求和指标。位置信息沿用原生体素坐标，Angle 来自原生导出的弧度值；脚本未做毫米换算。

这是原始配对数据，尚未生成 Path Transformer 训练用的弧长对齐、局部标架偏移 u/v、截面角 theta，也未生成 SDF 特征。正式划分训练/验证/测试集时，应按场景或布局家族分组，不能把同一地图的候选拆到不同集合中。

## 6. 重新编译

源码更新后必须重新编译；包里的 exe 对应本次测试时的源码，不会自动跟随原目录更新。

使用 Python 3.9+、CMake 3.16+ 和可用的 C++17/OpenMP 工具链。以下在已配置 MSVC 与 Ninja 的 Developer PowerShell 中运行：

```powershell
python .\build.py --busbar-source "E:\Learning\Subjects\路径生成\busbar_hybrid_state\busbar_hybrid_state" --voxel-source "E:\Learning\Subjects\路径生成\VP\Voxel_Planner" --build-dir "C:\Users\ASUS\Documents\Codex\busbar_test_build" --generator Ninja
```

新程序位于指定构建目录的 `bin/hybrid_validation_cli.exe`，随后通过 `run_tests.py --exe` 使用它。

构建目录必须是新目录，建议采用英文路径。构建器复制参与构建的源文件到该目录，并记录逐文件 SHA-256，不修改 E 盘原工程。两个项目原来的同名 demo 目标不直接混合；本 CMake 只额外构建师姐原始主程序及其库依赖。

本次测试遇到的环境问题：旧 MinGW 8.1 对中文路径及 std::filesystem 支持存在实际编译问题；MSVC 默认旧 OpenMP 模式不接受库中的无符号循环变量，因此本 CMake 选择 `/openmp:llvm`。当前 exe 已在这台机器运行验证。

若构建失败，查看构建目录的 `configure.log` 和 `build.log`。不要把编译或 DLL 问题记为算法无解。

## 7. 回归测试

测试 Python 结果处理逻辑：

```powershell
python -m unittest discover -s . -p test_runner.py -v
```

编译并运行原库的两项现有测试（把目录替换为自己的构建目录）：

```powershell
cmake --build "C:\Users\ASUS\Documents\Codex\busbar_test_build" --target test_sdf test_elastic_band --config Release
ctest --test-dir "C:\Users\ASUS\Documents\Codex\busbar_test_build" -C Release --output-on-failure
```

集成主程序没有调用 `SqpProjector`；其相关行为由现有 `test_elastic_band` 单独检查。不能把本次集成结果解释为已经验证原 `generate_dataset` 的全部产出。

## 8. 已发现但未改动的原项目行为

1. **端部障碍豁免**：Voxel_Planner 中有端部豁免机制。起点被占据的固定案例返回了粗路径，但师姐细规划器拒绝了它。因此脚本不会把粗规划成功单独作为入库依据。
2. **端部法向未传递**：师姐 `src/main.cpp` 构造 `TerminalConstraints` 时设置了切向与端部长度，没有赋值 start_normal/end_normal，粗规划使用 API 默认法向。本次 15 个样例均采用默认 Z 法向；其他法向需专门核实，不能根据本次结果声称覆盖。
3. **旧 path_validator.py 的含义**：它主要比较两个程序的状态、代价、节点数与路径距离，甚至两边都失败时也可能输出 PASS。它不是独立的全尺寸制造可行性证明，本入口没有把其 PASS 当作样本入库依据。
4. **工艺规则是否启用**：需在清单里提供企业确认的折弯 case 参数，并确认相关段类型开关。当前合成案例没有企业标准参数。

下一步需确认：你们所说“标准解”究竟以哪个检查接口及约束集合为准，以及真实地图、端点、尺寸、工艺参数和数据集目标格式。确认后可沿用本入口补充真实案例，而无需重新搭建调用链。
