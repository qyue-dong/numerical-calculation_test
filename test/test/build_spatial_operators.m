function ops = build_spatial_operators(N, h)
%% build_spatial_operators.m (Final Optimized Version)
% 针对 Dirichlet 边界条件优化的 4 阶紧致算子生成器
% 核心改进：使用“幽灵点法”处理近边界点，消除 2D 角点误差

I_N = speye(N);
e   = ones(N,1);

%% 1. 一维 4 阶一阶导 C1 (用于对流项)
% f'(x_i) ≈ (-f_{i+2} + 8f_{i+1} - 8f_{i-1} + f_{i-2}) / 12h
d0 = zeros(N,1);
d1 = (2/3)*ones(N,1);
d2 = (-1/12)*ones(N,1);
C1 = spdiags([-d2, -d1, d0, d1, d2], -2:2, N, N);

% 边界修正 (显式 4 阶偏心格式)
C1(1,1:5) = [-25, 48, -36, 16, -3]/12;
C1(2,1:5) = [-3, -10, 18, -6, 1]/12;
C1(N,N-4:N) = -fliplr(C1(1,1:5));
C1(N-1,N-4:N) = -fliplr(C1(2,1:5));

C1 = C1 / h;

%% 2. 一维 4 阶二阶导 C2_D (Dirichlet 专用 - 核心修正)
% 核心逻辑：利用 u_0 = -u_2 (反对称幽灵点) 消除边界误差
d0 = (-5/2)*ones(N,1);
d1 = (4/3)*ones(N,1);
d2 = (-1/12)*ones(N,1);
C2_D = spdiags([d2, d1, d0, d1, d2], -2:2, N, N);

% --- 关键修正开始 ---
% 左边界附近 (i=2):
% 标准中心格式: (-u_0 + 16u_1 - 30u_2 + 16u_3 - u_4) / 12h^2
% 代入 u_1=0 (Dirichlet), u_0 = -u_2 (幽灵点)
% 结果: (-(-u_2) + 0 - 30u_2 + 16u_3 - u_4) -> (-29u_2 + 16u_3 - u_4)
% 对应系数: [0, -29, 16, -1, 0]/12
C2_D(2, 1:5) = [16, -29, 16, -1, 0] / 12;

% 右边界附近 (i=N-1): 对称处理
C2_D(N-1, N-4:N) = fliplr(C2_D(2, 1:5));

% 边界点 (i=1, i=N): 
% 这些行会被 apply_bc_matrix_lhs 覆盖，但为了矩阵性质设为单位行
C2_D(1,:) = 0; C2_D(1,1) = 1; 
C2_D(N,:) = 0; C2_D(N,N) = 1; 
% --- 关键修正结束 ---

C2_D = C2_D / h^2;

%% 3. 一维 2 阶二阶导 C2_N (辅助用)
e2 = ones(N,1);
C2_N = spdiags([e2, -2*e2, e2], -1:1, N, N);
C2_N(1,1:2) = [-1, 1]; C2_N(N,N-1:N) = [1, -1];
C2_N = C2_N / h^2;

%% 4. 紧致算子系数矩阵 (B_i)
% 对应论文系数：1/12 * [1, 10, 1]
C3 = spdiags([e, 10*e, e], -1:1, N, N);
C4 = spdiags([-e, 14*e, -e], -1:1, N, N);
C5 = spdiags([-e, 8*e, -e], -1:1, N, N);

%% 5. 组装 2D 算子 (Kronecker Product)
ops.A1 = kron(C1, I_N);    % d/dx
ops.A2 = kron(I_N, C1);    % d/dy

ops.A3 = kron(C2_D, I_N);  % d2/dx2
ops.A5 = kron(I_N, C2_D);  % d2/dy2

% A4: 修正项 (h^2/12 * Delta)
ops.A4 = (h^2/12) * (ops.A3 + ops.A5);

% B_i 算子
ops.B2 = kron(C3, I_N) / 12;
ops.B3 = kron(I_N, C3) / 12;
ops.B4 = kron(C4, I_N) / 12;
ops.B5 = kron(I_N, C5) / 6;

%% 6. 组合 NS-alpha 核心矩阵
% H1 = B2 + B3 + A4
ops.H1 = (ops.B2 + ops.B3) + ops.A4;

% H2 = B2*A3 + B3*A5
ops.H2 = ops.B2 * ops.A3 + ops.B3 * ops.A5;

% 压力梯度
ops.Grad_x = ops.B4 * ops.A1;
ops.Grad_y = ops.B5 * ops.A2;

% 辅助拉普拉斯
ops.Lp_N = kron(C2_N, I_N) + kron(I_N, C2_N);

% 保存网格信息
ops.N = N;
ops.h = h;
end
