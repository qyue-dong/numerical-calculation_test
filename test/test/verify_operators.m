%% 1. 设置测试环境 (Setup)
clear; clc; close all;

% --- 定义测试参数 ---
N = 41;         % 网格点数
L = 1.0;        % 物理域长度
h = L / (N-1);  % 空间步长
x = linspace(0, L, N)'; % 创建一维网格列向量

% --- 定义测试函数及其解析导数 ---
f_analytic        = sin(2*pi*x);
df_dx_analytic    = 2*pi*cos(2*pi*x);
d2f_dx2_analytic  = -(2*pi)^2 * sin(2*pi*x);

fprintf('--- Verification Started for N=%d ---\n', N);

%% 2. 执行与计算 (Execute & Compute)

% --- 调用待验证的函数 ---
alpha_test = h; 
try
    ops = build_spatial_operators(N, h, alpha_test);
    fprintf('`build_spatial_operators` executed successfully.\n\n');
catch ME
    error('Failed to execute `build_spatial_operators`. Error: %s', ME.message);
end

% --- 临时的解决方案：直接在验证脚本里重新构建1D算子进行测试 ---
e  = ones(N, 1);
% C1 for test
C1_test = spdiags([-e, e], [-1, 1], N, N);
C1_test(1, 1:3) = [-3, 4, -1];
C1_test(N, N-2:N) = [1, -4, 3];
C1_test = C1_test / (2*h);

% C2 for test (Corrected block)
C2_test = spdiags([e, -2*e, e], [-1, 0, 1], N, N);
C2_test(1, 1:4) = [2, -5, 4, -1];
C2_test(N, N-3:N) = [-1, 4, -5, 2];
C2_test = C2_test / h^2;

% --- 计算数值导数 ---
df_dx_numeric = C1_test * f_analytic;
d2f_dx2_numeric = C2_test * f_analytic;


%% 3. 对比与分析 (Compare & Analyze)
% --- 验证 C1 (一阶导数) ---
figure('Name', 'Verification of C1 (First Derivative)');
subplot(1, 2, 1);
plot(x, df_dx_analytic, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytic Derivative');
hold on;
plot(x, df_dx_numeric, 'ro--', 'LineWidth', 1, 'DisplayName', 'Numeric Derivative (C1)');
grid on;
title('C1: Analytic vs. Numeric First Derivative');
xlabel('x'); ylabel('f''(x)'); legend('show');
subplot(1, 2, 2);
plot(x, abs(df_dx_numeric - df_dx_analytic), 'k-p', 'LineWidth', 1);
grid on;
title('Absolute Error of C1');
xlabel('x'); ylabel('|Error|');
error_C1 = max(abs(df_dx_numeric - df_dx_analytic));
fprintf('Max error for C1 (1st derivative): %e\n', error_C1);

% --- 验证 C2 (二阶导数) ---
figure('Name', 'Verification of C2 (Second Derivative)');
subplot(1, 2, 1);
plot(x, d2f_dx2_analytic, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytic Derivative');
hold on;
plot(x, d2f_dx2_numeric, 'ro--', 'LineWidth', 1, 'DisplayName', 'Numeric Derivative (C2)');
grid on;
title('C2: Analytic vs. Numeric Second Derivative');
xlabel('x'); ylabel('f''''(x)'); legend('show');
subplot(1, 2, 2);
plot(x, abs(d2f_dx2_numeric - d2f_dx2_analytic), 'k-p', 'LineWidth', 1);
grid on;
title('Absolute Error of C2');
xlabel('x'); ylabel('|Error|');
error_C2 = max(abs(d2f_dx2_numeric - d2f_dx2_analytic));
fprintf('Max error for C2 (2nd derivative): %e\n', error_C2);

fprintf('\n--- Verification Finished ---\n');

% --- 可选：检查通过 build_spatial_operators 生成的2D算子稀疏结构 ---
figure('Name', 'Sparsity Patterns of 2D Operators');
subplot(2, 3, 1); spy(ops.A1); title('Sparsity of A1 (d/dx)');
subplot(2, 3, 2); spy(ops.A2); title('Sparsity of A2 (d/dy)');
subplot(2, 3, 3); spy(ops.Lp); title('Sparsity of Laplacian');
subplot(2, 3, 4); spy(ops.H1); title('Sparsity of H1 (Mass)');
subplot(2, 3, 5); spy(ops.H2); title('Sparsity of H2 (Viscosity)');
subplot(2, 3, 6); spy(ops.M_filter); title('Sparsity of M_{filter}');