function ops = build_spatial_operators(N, h, alpha)

fprintf('   - Building spatial operators for N=%d grid...\n', N);

%% 1. 构建基础的一维(1D)微分矩阵 (1D Primitives for Dirichlet BCs)
I1 = speye(N);
e  = ones(N, 1);

% --- C1: 一阶导数算子 (First Derivative Operator) ---
% 内部使用二阶中心差分，边界使用二阶单边差分
C1 = spdiags([-e, e], [-1, 1], N, N);
C1(1, 1:3) = [-3, 4, -1];
C1(N, N-2:N) = [1, -4, 3];
C1 = C1 / (2*h);

% --- C2: 二阶导数算子 (Second Derivative Operator) ---
% 内部使用二阶中心差分，边界使用二阶单边差分
C2 = spdiags([e, -2*e, e], [-1, 0, 1], N, N);
% 在边界处使用二阶精度的单边差分
C2(1, 1:4) = [2, -5, 4, -1]; 
C2(N, N-3:N) = [-1, 4, -5, 2]; % <--- 系数[-1, 4, -5, 2]是正确的
C2 = C2 / h^2;

% --- 【新步骤】为B算子创建一个内部版本的C2 ---
% 这个版本的C2在边界处是0，确保B算子在边界处是单位算子
C2_for_B = spdiags([e, -2*e, e], [-1, 0, 1], N, N);
C2_for_B(1, :) = 0;   % 第一行全设为0
C2_for_B(N, :) = 0;   % 最后一行全设为0
C2_for_B = C2_for_B / h^2;

% --- B_i 算子的一维构建块 ---
% 【修正】使用 C2_for_B 来构建 B 算子
B1_1D = I1 + (h^2/12) * C2_for_B; 
B2_1D = (1/12) * spdiags([e, 10*e, e], [-1, 0, 1], N, N);
B2_1D(1,1) = 1; B2_1D(1,2)=0;
B2_1D(N,N) = 1; B2_1D(N,N-1)=0;
B3_1D = B2_1D;
% 【修正】使用 C2_for_B 来构建 B 算子
B4_1D = I1 - (h^2/6) * C2_for_B;
B5_1D = B4_1D;

%% 2. 使用Kronecker积构建2D算子 (Construct 2D Operators via Kronecker Product)
% 修正: kron(A, I) -> x方向; kron(I, A) -> y方向

% --- A_i 算子 (A_i Operators) ---
ops.A1 = kron(C1, I1);       % A1 = d/dx
ops.A2 = kron(I1, C1);       % A2 = d/dy
ops.A3 = kron(C2, I1);       % A3 = d^2/dx^2
ops.A5 = kron(I1, C2);       % A5 = d^2/dy^2
ops.Lp = ops.A3 + ops.A5;    % 标准拉普拉斯算子

% --- H_i 算子，用于投影法 (H_i Operators for Projection Method) ---
ops.H1 = kron(B1_1D, B1_1D); % 质量矩阵
ops.H2 = kron(B2_1D, I1) * ops.A3 + kron(I1, B3_1D) * ops.A5; % 粘性矩阵

%% 3. 组合构建核心算子 (Core Combined Operators)

% --- 滤波速度更新算子 (Filter Update Operator) ---
% 修正: 'alpha' 来自函数输入
ops.M_filter = ops.H1 - alpha^2 * ops.H2;

% --- 压力相关的算子 (Pressure-related Operators) ---
ops.Grad_x = kron(B4_1D, I1) * ops.A1;
ops.Grad_y = kron(I1, B5_1D) * ops.A2;
ops.Div = ops.Grad_x + ops.Grad_y; % 高阶散度算子

fprintf('   - All spatial operators constructed successfully.\n\n');

end