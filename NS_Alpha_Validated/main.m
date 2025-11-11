%% ==========================================================
%% >>>>> 主程序: v1.0 - 可靠的二阶收敛测试 <<<<<
%% ==========================================================
clear; clc; close all;

%% 1. 参数设置
nu = 0.1;
alpha = 0.1; % 尽管在此线性测试中alpha不起作用
T_final = 0.1;

refinement_levels = 4;
N_base = 8;
N_list = N_base * 2.^(0:refinement_levels-1); % -> [8, 16, 32, 64]
errors_L2 = zeros(1, refinement_levels);

fprintf('=== 开始【v1.0】收敛阶测试 (纯线性、二阶) ===\n');

%% 2. 主循环
for i = 1:length(N_list)
    N = N_list(i);
    h = (2*pi) / N;
    fprintf('\n--- 正在处理网格: N = %d, h = %f ---\n', N, h);

    % --- A. 构建网格和算子 ---
    ops = build_operators(N, h);
    vec_p = linspace(0, 2*pi, N+1); vec_p = vec_p(1:end-1);
    [X, Y] = meshgrid(vec_p);
    
    % --- B. 时间步长 ---
    tau = 0.1 * h^2;
    num_steps = ceil(T_final / tau);
    tau = T_final / num_steps;

    % --- C. 初始化 ---
    mms = ManufacturedSolution(nu, alpha);
    [V, W] = mms.getInitialCondition(X, Y); % t=0

    % --- D. 时间推进循环 ---
    fprintf('时间步长 tau = %e, 总步数 = %d\n', tau, num_steps);
    for n = 1:num_steps
        t_center = (n - 0.5) * tau;
        
        % 获取源项
    
    [Sv_full, Sw_full] = mms.getSourceTerms(t_center, X, Y);
    
    % 为了求解线性的斯托克斯方程，我们必须从源项中减去精确解的非线性部分
    [V_center, W_center, ~, V_bar_center, W_bar_center] = mms.getSolution(t_center, X, Y);
    [NL_v_center, NL_w_center] = compute_nonlinear_term(V_center, W_center, V_bar_center, W_bar_center, ops);

    Sv_matrix = Sv_full - reshape(NL_v_center, N, N);
    Sw_matrix = Sw_full - reshape(NL_w_center, N, N);
        

        [Sv, Sw] = mms.getSourceTerms(t_center, X, Y);

        % 调用求解器
        [V, W] = Your_NS_Alpha_Solver(V, W, h, tau, nu, Sv, Sw, ops);
    end

    % --- E. 计算误差 ---
    [V_exact, W_exact] = mms.getSolution(T_final, X, Y);
    error_v = V - V_exact;
    error_w = W - W_exact;
    errors_L2(i) = sqrt(h^2 * sum(error_v(:).^2 + error_w(:).^2));
    
    fprintf('计算完成。L2误差 = %e\n', errors_L2(i));
end

%% 3. 计算并展示收敛阶
orders_L2 = log(errors_L2(1:end-1) ./ errors_L2(2:end)) / log(2);
fprintf('\n==================== 最终收敛结果 ====================\n');
fprintf('%-10s %-15s %-15s\n', 'N', 'L2 Error', 'L2 Order');
fprintf('%-10d %-15.4e %-15s\n', N_list(1), errors_L2(1), '---');
for i = 2:length(orders_L2)
    fprintf('%-10d %-15.4e %-15.4f\n', N_list(i), errors_L2(i), orders_L2(i-1));
end
fprintf('%-10d %-15.4e %-15.4f\n', N_list(end), errors_L2(end), orders_L2(end));