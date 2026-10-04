# Path Transformer 最小推理部署

此目录提供 Voxel_Planner 粗路径的 Path Transformer 推理接口。输入为集成测试产生的 `accepted.jsonl`；每个路径点提供体素坐标和法向、切向。模型输出局部法向偏移 `u`、副法向偏移 `v` 和截面旋转角 `theta`（弧度），并重建中心线。端点修正固定为零。

当前项目没有训练权重。`--smoke` 仅运行未训练模型，输出恒等路径以验证输入、模型前向和文件写出，状态为 `UNTRAINED_SMOKE_ONLY`；这些路径不是模型优化结果。真实推理必须提供包含 `architecture: PathTransformer-v1`、`trained: true` 和 `state_dict` 的兼容 checkpoint，输出状态仍为 `CHECKPOINT_INFERENCE_UNVALIDATED`，直到完成约束校验。推理结果不自动交给制造流程。

使用 Python 3.9 和 PyTorch CPU 版即可运行；在项目根目录执行：

```powershell
python -m pip install numpy==1.26.4
python -m pip install torch==2.2.2 --index-url https://download.pytorch.org/whl/cpu
python .\path_transformer\infer.py --input .\weekly_2026-09-27\verification_batch_2026-09-27\accepted.jsonl --output .\weekly_2026-09-27\my_smoke.jsonl --smoke
```

输出记录 `predicted_centerline`、`local_corrections_u_v_theta`、`hard_constraints_verified: false`。输入特征依次为相对路径位置、切向、法向、归一化弧长、宽度、厚度和间隙占位值，共 13 维。现有集成导出没有逐点 SDF，间隙特征暂填 0，并在输出标注 `clearance_feature_available: false`。当前使用 2 层、4 头、64 维的 PyTorch Transformer Encoder；这是用于验证部署接口的最小结构，尚无训练实验或路径质量结论。
