function ops = build_spatial_operators(N, h, alpha)

%% 1. 构建基础的一维(1D)微分矩阵 (Construct 1D Primitives)
% -------------------------------------------------------------------------
% 这些是构建2D算子的基础模块。我们使用 spdiags 来高效创建稀疏矩阵。
% -------------------------------------------------------------------------

I1 = speye(N);      % N x N 的稀疏单位矩阵
e  = ones(N, 1);    % N x 1 的全1向量

% --- C1: 对应论文中 A1, A2 的核心，标准二阶精度一阶中心差分 ---
% d/dx f_i ≈ (f_{i+1} - f_{i-1}) / (2h)

C1 = spdiags([-e, e], [-1, 1], N, N);

% 周期性边界条件 (Periodic Boundary Conditions)
C1(1, N) = -1; % f(1)的导数需要 f(N)
C1(N, 1) =  1; % f(N)的导数需要 f(1)
C1 = C1 / (2*h); % 整个矩阵除以 2h

% --- C2: 对应论文中 A3, A5 的核心，标准二阶精度二阶中心差分 ---
% d^2/dx^2 f_i ≈ (f_{i+1} - 2f_i + f_{i-1}) / h^2

C2 = spdiags([e, -2*e, e], [-1, 0, 1], N, N);

% 周期性边界条件 (Periodic Boundary Conditions)
C2(1, N) = 1; % f(1)的二阶导数需要 f(N)
C2(N, 1) = 1; % f(N)的二阶导数需要 f(1)
C2 = C2 / h^2;

% 对于Dirichlet边界条件，这个标准格式是足够的。

% --- B_i 算子的一维构建块 (1D building blocks for B_i operators) ---
% B1_1D: 对应算子 (I + (h^2/12)*d^2/dx^2)

B1_1D = I1 + (h^2/12) * C2;

% B2_1D / B3_1D: 对应论文中 (1/12)[1, 10, 1] 的紧致算子

B2_1D = (1/12) * spdiags([e, 10*e, e], [-1, 0, 1], N, N);

% 边界修正，使其在边界上退化为单位算子

B2_1D(1,1) = 1; B2_1D(1,2)=0;
B2_1D(N,N) = 1; B2_1D(N,N-1)=0;
B3_1D = B2_1D; % 在这个模型中 B2 和 B3 相同

% B4_1D / B5_1D: 对应算子 (I - (h^2/6)*d^2/dx^2)，这是四阶梯度/散度的修正项
B4_1D = I1 - (h^2/6) * C2;
B5_1D = B4_1D; % B4 和 B5 相同

%% 2. 使用Kronecker积构建2D算子 (Construct 2D Operators)
% -------------------------------------------------------------------------
% kron(A, I) 将算子A应用于网格的每一列 (x方向)
% kron(I, A) 将算子A应用于网格的每一行 (y方向)
% -------------------------------------------------------------------------

% --- A_i 算子 (A_i Operators) ---
ops.A1 = kron(I1,C1);   % A1 = d/dx
ops.A2 = kron(C1,I1);   % A2 = d/dy
ops.A3 = kron(I1,C2);   % A3 = d^2/dx^2
ops.A5 = kron(C2,I1);   % A5 = d^2/dy^2
ops.Lp = ops.A3 + ops.A5; % 标准拉普拉斯算子 (Standard Laplacian)

% --- H_i 算子，用于投影法 (H_i Operators for Projection Method) ---
% H1: 质量矩阵 (Mass Matrix), 来自 B1 算子的2D形式
ops.H1 = kron(B1_1D, B1_1D);

% H2: 粘性矩阵 (Viscosity Matrix), 来自 B2, B3 和 A3, A5
ops.H2 = kron(B2_1D, I1) * ops.A3 + kron(I1, B3_1D) * ops.A5;

%% 3. 组合构建投影法所需的核心算子 (Core Operators for Projection Method)
% -------------------------------------------------------------------------
% 这些是时间推进循环中将要直接使用的组合算子。
% -------------------------------------------------------------------------

% --- 滤波速度更新算子 (Filter Update Operator) ---
% M_filter = H1 - alpha^2 * H2

ops.M_filter = ops.H1 - alpha^2 * ops.H2;

% --- 压力相关的算子 (Pressure-related Operators) ---
% Grad_x, Grad_y: 高阶梯度算子 (High-order Gradient)

ops.Grad_x = kron(B4_1D, I1) * ops.A1;
ops.Grad_y = kron(I1, B5_1D) * ops.A2;

% Div: 高阶散度算子 (High-order Divergence)
% 注意: 散度是梯度的负共轭，这里我们直接构建用于压力泊松方程的算子

ops.Div = ops.Grad_x + ops.Grad_y;

fprintf('   - All spatial operators constructed successfully.\n\n');

end