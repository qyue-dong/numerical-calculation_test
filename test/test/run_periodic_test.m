function run_periodic_test()
%% 1D-like 2D Poisson Test (Final Fixed Version)
% 目的：消除角点干扰 + 修复网格步长定义的冲突
% 修正点：X方向 hx = L/N (周期); Y方向 hy = L/(N-1) (Dirichlet)

clear; clc; close all;

L = 1.0;
grid_sizes = [16, 24, 32, 48, 64]; 
errors = zeros(length(grid_sizes), 1);

fprintf('============================================================\n');
fprintf('Periodic-Dirichlet Poisson Verification (h_x != h_y)\n');
fprintf('============================================================\n');

for k = 1:length(grid_sizes)
    N = grid_sizes(k);
    
    % --- 关键修正：分离步长 ---
    hx = L / N;       % 周期性：N 个点不包含终点 L
    hy = L / (N - 1); % Dirichlet：N 个点包含终点 L
    
    % --- 1. 构建算子 ---
    I = speye(N);
    e = ones(N,1);
    
    % A. X方向：周期性 4 阶二阶导 (Dxx_Per)
    % 系数: [-1/12, 4/3, -5/2, 4/3, -1/12]
    d0 = (-5/2)*e; d1 = (4/3)*e; d2 = (-1/12)*e;
    Dxx_Per = spdiags([d2 d1 d0 d1 d2], -2:2, N, N);
    % 补角 (周期性 wrap-around)
    Dxx_Per(1, N-1:N) = [-1/12, 4/3]; Dxx_Per(2, N) = -1/12;
    Dxx_Per(N, 1:2) = [4/3, -1/12];   Dxx_Per(N-1, 1) = -1/12;
    
    Dxx_Per = Dxx_Per / hx^2; % <--- 使用 hx
    
    % B. Y方向：Dirichlet 4 阶二阶导 (Dyy_Dir)
    % 系数同上，但边界不同
    C2_D = spdiags([d2, d1, d0, d1, d2], -2:2, N, N);
    
    % Ghost Point 修正 (u0 = -u2) -> 适用于 sin(pi*y) 这种边界为0且反对称的函数
    % 系数: [16, -29, 16, -1, 0] / 12
    C2_D(2, 1:5) = [16, -29, 16, -1, 0] / 12; 
    C2_D(N-1, N-4:N) = fliplr(C2_D(2, 1:5));
    
    % 边界行置位 (稍后会被 apply_bc_lifting 覆盖，但保持矩阵结构完整)
    C2_D(1,:) = 0; C2_D(1,1) = 1; 
    C2_D(N,:) = 0; C2_D(N,N) = 1; 
    
    Dyy_Dir = C2_D / hy^2; % <--- 使用 hy
    
    % C. 组装 2D 拉普拉斯
    % X是Outer，Y是Inner
    % Lp = Dxx(X) + Dyy(Y)
    Lp = kron(Dxx_Per, I) + kron(I, Dyy_Dir);
    
    % --- 2. 构造真解 ---
    % X: 0 ~ L-hx (周期性，不含 L)
    % Y: 0 ~ L    (Dirichlet，含 L)
    x_1d = linspace(0, L-hx, N); 
    y_1d = linspace(0, L, N);   
    [X, Y] = meshgrid(x_1d, y_1d);
    
    % 真解: u = sin(2*pi*x) * sin(pi*y)
    u_ex = sin(2*pi*X) .* sin(pi*Y);
    
    % 源项 f = -Delta u = -(u_xx + u_yy)
    % u_xx = -4*pi^2 * u
    % u_yy = -pi^2 * u
    % f = -(-4pi^2 - pi^2)u = 5*pi^2 * u
    f_ex = 5 * pi^2 * u_ex;
    
    % --- 3. 边界处理 ---
    % Dirichlet 边界在 Y=0 (Row 1) 和 Y=L (Row N)
    idx_bot = 1:N:(N*N);
    idx_top = N:N:(N*N);
    idx_bd = [idx_bot, idx_top]';
    
    % 应用边界 (Lifting)
    u_bd_vals = u_ex(idx_bd); % 理论上全是 0
    
    A_sys = -Lp;
    b_rhs = f_ex(:);
    
    [A_BC, b_mod] = apply_bc_lifting(A_sys, b_rhs, idx_bd, u_bd_vals);
    
    % --- 4. 求解 ---
    u_sol = A_BC \ b_mod;
    
    % --- 5. 误差 ---
    err = norm(u_sol - u_ex(:)) / norm(u_ex(:));
    errors(k) = err;
    
    fprintf('Grid %dx%d: Error = %.4e\n', N, N, err);
end

% --- 绘图 ---
h_vals = L ./ grid_sizes; % 作图用的参考 h
figure('Color','w');
loglog(h_vals, errors, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog(h_vals, errors(1)*(h_vals/h_vals(1)).^4, 'k--', 'LineWidth', 1.5, 'DisplayName', '4th Order Ref');
grid on; xlabel('h (approx)'); ylabel('L2 Relative Error'); 
title('Periodic-Dirichlet 4th Order Verification');
legend('Our Code', '4th Order Slope', 'Location', 'best');

grid_ratios = grid_sizes(2:end) ./ grid_sizes(1:end-1);
error_ratios = errors(1:end-1) ./ errors(2:end);
% Rate = log(Error_Ratio) / log(Grid_Ratio)
rates = log(error_ratios) ./ log(grid_ratios');

fprintf('\n--- Corrected Convergence Rates ---\n');
disp(rates);

end

function [A_mod, b_mod] = apply_bc_lifting(A, b, idx_bd, u_bd_vals)
    [n, ~] = size(A);
    u_known = zeros(n, 1);
    u_known(idx_bd) = u_bd_vals;
    
    % 移项修正：把边界值对内部的贡献移到 RHS
    b_mod = b - A * u_known;
    
    I = speye(n);
    d = ones(n, 1); d(idx_bd) = 0; 
    D = spdiags(d, 0, n, n);
    
    % 矩阵置1
    A_mod = D * A + (I - D);
    
    % RHS 边界强制赋值
    b_mod(idx_bd) = u_bd_vals;
end