%% 1. 初始化环境 (Initialize Environment)
clear;         % 清除工作区中的所有变量
clc;           % 清理命令窗口
close all;     % 关闭所有图形窗口

fprintf('--- NS-alpha Model Solver: Framework Initialization ---\n');
fprintf('Stage: Day 1 (Morning) - Setting up the simulation environment.\n\n');


%% 2. 参数定义区 (Parameter Definition Section)
% -------------------------------------------------------------------------
% 将所有可调参数集中放置在此处，方便管理、修改和复现实验。
% -------------------------------------------------------------------------
fprintf('2. Defining simulation parameters...\n');

% --- 物理参数 (Physical Parameters) ---
Re = 1000;          % 雷诺数 (Reynolds Number, Re)
L  = 1.0;           % 方腔边长 (Domain Length)
nu = 1 / Re;        % 运动粘度 (Kinematic Viscosity, ν), 根据Re自动计算

% --- 数值参数 (Numerical Parameters) ---
% 这是为最终生成 "参考解" 设置的参数
N  = 256;           % 单方向的网格点数 (Number of grid points in one direction)
h  = L / N;         % 空间步长 (Spatial step size, h), 根据N和L自动计算
dt = 1e-5;          % 时间步长 (Time step, τ), 为参考解选择一个非常小的值
T_final = 100;      % 模拟总时间 (Total simulation time)

% --- NS-alpha 模型特定参数 (Model-Specific Parameters) ---
alpha = h;          % 滤波半径 (Filter radius, α), ตามแผนงาน设定为与网格尺寸关联

% --- 控制与求解器参数 (Control & Solver Parameters) ---
tolerance = 1e-10;  % 稳态判定的残差容忍度 (Residual tolerance for steady-state)
plot_interval = 5000;% 每隔多少个时间步更新一次可视化结果 (Visualization update interval)

% --- 在命令窗口打印关键参数以供核对 ---
fprintf('   - Reynolds Number (Re): %d\n', Re);
fprintf('   - Grid Resolution: %d x %d\n', N, N);
fprintf('   - Spatial Step (h): %.4f\n', h);
fprintf('   - Time Step (dt): %.1E\n', dt);
fprintf('   - Filter Radius (alpha): %.4f\n\n', alpha);


%% 3. 网格生成 (Mesh Generation)
% -------------------------------------------------------------------------
% 创建计算所需的二维笛卡尔坐标网格。
% 我们使用单元中心网格(cell-centered grid)，即网格点位于单元的中心。
% -------------------------------------------------------------------------
fprintf('3. Generating computational mesh...\n');

% 定义每个方向上单元中心的坐标
x_coords = linspace(h/2, L-h/2, N);
y_coords = linspace(h/2, L-h/2, N);

% 使用 meshgrid 创建二维坐标矩阵
[X, Y] = meshgrid(x_coords, y_coords);

fprintf('   - Mesh generated successfully.\n\n');


%% 4. 场变量初始化 (Field Variable Initialization)
% -------------------------------------------------------------------------
% 初始化所有需要求解的物理场变量。根据顶盖驱动方腔流的初始条件，
% 所有场在 t=0 时刻均为零。
% -------------------------------------------------------------------------
fprintf('4. Initializing field variables...\n');

% v: 速度场 u 的第一个分量 (x-velocity component)
v = zeros(N, N);

% w: 速度场 u 的第二个分量 (y-velocity component)
w = zeros(N, N);

% v_bar: 滤波后速度场 u_bar 的第一个分量
v_bar = zeros(N, N);

% w_bar: 滤波后速度场 u_bar 的第二个分量
w_bar = zeros(N, N);

% p: 压力场 (Pressure field)
p = zeros(N, N);

fprintf('   - All fields (v, w, v_bar, w_bar, p) initialized to zero.\n\n');

%% 5. 构建空间离散算子 (Construct Spatial Operators)
% -------------------------------------------------------------------------
% 调用独立的函数来构建所有的空间算子矩阵。
% 所有算子都将存储在名为 'ops' 的结构体中。
% -------------------------------------------------------------------------
ops = build_spatial_operators(N, h, alpha);
