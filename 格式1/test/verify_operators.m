
%验证build_spatial_operators.m中矩阵算子构建的正确性

%% 1. 设置测试环境
clear; clc; close all;

% 定义一系列网格尺寸用于收敛性分析
% 建议使用 N = 2^k + 1 的形式，这样 h 减半时网格点能重合，且计算方便
N_list = [21, 41, 81, 161]; 

error_1st_derivative = zeros(length(N_list), 1);
error_2nd_derivative = zeros(length(N_list), 1);
h_list = zeros(length(N_list), 1);

fprintf('========= 开始四阶精度验证 (基于提供的显式算子 A1, A3) =========\n');

for k = 1:length(N_list)
    N = N_list(k);
    L = 1.0;
    h = L / (N-1); % 确保 h 定义与您函数内部一致
    h_list(k) = h;
    
    % --- 生成 2D 网格 ---
    % 对应 meshgrid: x 沿列变化 (第二维), y 沿行变化 (第一维)
    x_1d = linspace(0, L, N);
    y_1d = linspace(0, L, N);
    [X, Y] = meshgrid(x_1d, y_1d);
    
    % 将网格拉直为列向量 (Stack columns)
    x = X(:);
    y = Y(:);
    
    % --- 定义测试函数 (2D) ---
    % 选用光滑周期函数，避免边界非齐次干扰，专测核心算子精度
    % f(x,y) = sin(2*pi*x) * sin(2*pi*y)
    f_exact = sin(2*pi*x) .* sin(2*pi*y);
    
    % --- 定义解析导数 (2D) ---
    % df/dx (对应 A1)
    df_dx_exact = 2*pi * cos(2*pi*x) .* sin(2*pi*y);
    
    % d^2f/dx^2 (对应 A3)
    d2f_dx2_exact = -(2*pi)^2 * sin(2*pi*x) .* sin(2*pi*y);
    
    % --- 调用您的算子构建函数 ---
    % 注意：您的函数定义为 build_spatial_operators(N, h)，不需要 alpha
    ops = build_spatial_operators(N, h);
    
    % --- 计算数值导数 ---
    % 根据您的代码：
    % ops.A1 = kron(C1, I_N)  -> 对应 d/dx (作用于 x 维度)
    % ops.A3 = kron(C2_D, I_N) -> 对应 d^2/dx^2 (作用于 x 维度)
    % 并且您的 C1, C2_D 是显式矩阵，因此直接相乘即可
    
    df_dx_num   = ops.A1 * f_exact;
    d2f_dx2_num = ops.A3 * f_exact;
    
    % --- 计算误差 (L_infinity Norm: 最大绝对误差) ---
    % 排除边界点影响（可选）：
    % 也就是只统计内部点的误差，因为边界点有时精度会降一阶，但整体仍收敛。
    % 这里我们先统计全局误差。
    error_1st_derivative(k) = norm(df_dx_num - df_dx_exact, inf);
    error_2nd_derivative(k) = norm(d2f_dx2_num - d2f_dx2_exact, inf);
    
    fprintf('Grid: %3d x %3d, h = %.5f | Err(d/dx)=%.4e | Err(d^2/dx^2)=%.4e\n', ...
            N, N, h, error_1st_derivative(k), error_2nd_derivative(k));
end

%% 2. 计算收敛阶 (Order of Accuracy)
% Order = log2( Error(Coarse) / Error(Fine) )
order_1st = log2(error_1st_derivative(1:end-1) ./ error_1st_derivative(2:end));
order_2nd = log2(error_2nd_derivative(1:end-1) ./ error_2nd_derivative(2:end));

fprintf('\n--- 收敛阶结果 (理论值应接近 4.0) ---\n');
fprintf('N的变化 \t\t 一阶导(A1)阶数 \t 二阶导(A3)阶数\n');
for k = 1:length(order_1st)
    fprintf('%3d -> %3d \t %.4f \t\t %.4f\n', ...
            N_list(k), N_list(k+1), order_1st(k), order_2nd(k));
end

%% 3. 绘图验证 (Log-Log Plot)
figure('Name', 'Spatial Operators Convergence Test', 'Color', 'w');

subplot(1,2,1);
loglog(h_list, error_1st_derivative, '-bo', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'A1 Error (d/dx)');
hold on;
% 添加斜率=4的参考线
ref_line = error_1st_derivative(end) * (h_list / h_list(end)).^4;
loglog(h_list, ref_line, 'k--', 'LineWidth', 1.5, 'DisplayName', '4th Order Slope');
grid on; axis tight;
xlabel('Grid spacing h'); ylabel('Max Error (L_\infty)');
title('1st Derivative (A1) Accuracy');
legend('Location', 'northwest');
set(gca, 'XDir', 'reverse'); % h 减小方向向右

subplot(1,2,2);
loglog(h_list, error_2nd_derivative, '-rs', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'A3 Error (d^2/dx^2)');
hold on;
% 添加斜率=4的参考线
ref_line_2 = error_2nd_derivative(end) * (h_list / h_list(end)).^4;
loglog(h_list, ref_line_2, 'k--', 'LineWidth', 1.5, 'DisplayName', '4th Order Slope');
grid on; axis tight;
xlabel('Grid spacing h'); ylabel('Max Error (L_\infty)');
title('2nd Derivative (A3) Accuracy');
legend('Location', 'northwest');
set(gca, 'XDir', 'reverse');