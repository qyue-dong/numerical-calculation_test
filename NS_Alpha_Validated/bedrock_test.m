%% =================================================================
%% >>>>>>>>>>>> 基岩测试 (Bedrock Test) —— 修复版 <<<<<<<<<<<<
%  一个完全独立的、无任何外部依赖的最终验证脚本
%  ✔️ 修复：不再手动构造源项，改用标准 fractional-step 方法
%  ✔️ 精确解：Taylor–Green 涡（NS 解析解）
%  ✔️ 时间推进：显式对流 + CN 扩散 + projection（一阶时间，二阶空间）
%  ✔️ 预期收敛：L2 误差 ~ O(h^2)（因 tau ~ h^2）
%% =================================================================
clear; clc; close all;

%% 1. 参数设置
nu = 0.1;
T_final = 0.1;
refinement_levels = 4;
N_base = 8;
N_list = N_base * 2.^(0:refinement_levels-1);
errors_L2 = zeros(1, refinement_levels);

fprintf('=== 开始【基岩测试】-- 修复版（标准 fractional-step） ===\n');

%% 2. 主循环
for i = 1:length(N_list)
    N = N_list(i);
    h = (2*pi) / N;
    fprintf('\n--- 正在处理网格: N = %d (h = %.3e) ---\n', N, h);
    N2 = N^2;

    % --- A. 构建周期性差分算子 ---
    I_N = speye(N);
    e = ones(N, 1);
    
    % 一阶导 (中心差分，周期性)
    D1 = spdiags([-e, e], [-1, 1], N, N); 
    D1(1, N) = -1; D1(N, 1) = 1;
    
    % 二阶导 (中心差分，周期性)
    D2 = spdiags([e, -2*e, e], [-1, 0, 1], N, N); 
    D2(1, N) = 1; D2(N, 1) = 1;
    
    % 拉普拉斯 & 梯度算子
    Dx = (1/(2*h)) * kron(I_N, D1);   % ∂/∂x
    Dy = (1/(2*h)) * kron(D1, I_N);   % ∂/∂y
    Laplacian = (1/h^2) * (kron(I_N, D2) + kron(D2, I_N)); % Δ

    % --- B. 时间步长设置（满足 parabolic CFL: tau ~ h^2）---
    tau = 0.1 * h^2;  % <= 0.5*h^2/nu for stability (nu=0.1 → 0.1*h^2 OK)
    num_steps = ceil(T_final / tau);
    tau = T_final / num_steps;  % 修正为整除

    % --- C. 网格与精确解定义 ---
    vec_p = linspace(0, 2*pi, N+1); vec_p = vec_p(1:end-1); % [0, 2π)
    [X, Y] = meshgrid(vec_p);
    
    % Taylor–Green 涡（NS 精确解）
    v_exact_func = @(t, x, y) -cos(x) .* sin(y) * exp(-2*nu*t);
    w_exact_func = @(t, x, y)  sin(x) .* cos(y) * exp(-2*nu*t);
    p_exact_func = @(t, x, y) -0.25 * (cos(2*x) + cos(2*y)) * exp(-4*nu*t);
    
    % 初始条件（t=0）
    V = v_exact_func(0, X, Y);
    W = w_exact_func(0, X, Y);

    % --- D. 时间推进循环（标准 fractional-step）---
    fprintf('时间步长 tau = %e, 总步数 = %d\n', tau, num_steps);
    
    for n = 1:num_steps
        % 当前步解向量化
        v_n = reshape(V, N2, 1);
        w_n = reshape(W, N2, 1);
        
        % 计算空间导数（用于对流项）
        Vx = reshape(Dx * v_n, N, N);
        Vy = reshape(Dy * v_n, N, N);
        Wx = reshape(Dx * w_n, N, N);
        Wy = reshape(Dy * w_n, N, N);
        
        % 网格场（用于点乘）
        V_grid = reshape(v_n, N, N);
        W_grid = reshape(w_n, N, N);
        
        % 显式计算非线性项 (u·∇)u at t = n*tau
        conv_v = V_grid .* Vx + W_grid .* Vy;   % (u·∇)v
        conv_w = V_grid .* Wx + W_grid .* Wy;   % (u·∇)w
        
        conv_v_vec = reshape(conv_v, N2, 1);
        conv_w_vec = reshape(conv_w, N2, 1);
        
        % --- 预测步：v* = v^n - τ*(u·∇v)^n + (τν/2)*Δ(v^n + v*) ---
        LHS = speye(N2) - 0.5 * tau * nu * Laplacian;
        RHS_v = v_n - tau * conv_v_vec + 0.5 * tau * nu * (Laplacian * v_n);
        RHS_w = w_n - tau * conv_w_vec + 0.5 * tau * nu * (Laplacian * w_n);
        
        v_star = LHS \ RHS_v;
        w_star = LHS \ RHS_w;
        
        % --- 压力求解：Δp = (1/τ) ∇·u* ---
        div_u_star = Dx * v_star + Dy * w_star;
        RHS_p = (1/tau) * div_u_star;
        RHS_p = RHS_p - mean(RHS_p);  % 相容性：∫RHS_p = 0
        
        p_next = Laplacian \ RHS_p;
        p_next = p_next - mean(p_next);  % 双重保险：p 零均值
        
        % --- 速度校正：u^{n+1} = u* - τ ∇p ---
        v_next = v_star - tau * (Dx * p_next);
        w_next = w_star - tau * (Dy * p_next);
        
        % 更新场
        V = reshape(v_next, N, N);
        W = reshape(w_next, N, N);
        
        % --- 【调试】监控散度（每 20% 步打印一次）---
        if mod(n, max(1, floor(num_steps/5))) == 0
            div_u = reshape(Dx * v_next + Dy * w_next, N, N);
            max_div = max(abs(div_u(:)));
            fprintf('  Step %3d/%d: max|∇·u| = %.2e\n', n, num_steps, max_div);
        end
    end

    % --- E. 计算最终 L2 误差 ---
    V_exact_tf = v_exact_func(T_final, X, Y);
    W_exact_tf = w_exact_func(T_final, X, Y);
    
    err_v = V - V_exact_tf;
    err_w = W - W_exact_tf;
    errors_L2(i) = sqrt(h^2 * sum(err_v(:).^2 + err_w(:).^2));
    
    fprintf('计算完成。L2误差 = %.4e\n', errors_L2(i));
end

%% 3. 收敛阶分析
fprintf('\n==================== 基岩测试结果（修复版） ====================\n');
fprintf('%-10s %-15s %-15s\n', 'N', 'L2 Error', 'L2 Order');
fprintf('%-10d %-15.4e %-15s\n', N_list(1), errors_L2(1), '---');
if length(N_list) > 1
    orders_L2 = log(errors_L2(1:end-1) ./ errors_L2(2:end)) / log(2);
    for i = 2:length(N_list)
        fprintf('%-10d %-15.4e %-15.4f\n', N_list(i), errors_L2(i), orders_L2(i-1));
    end
end

% 4. 【可选】绘图（若想可视化最后一层结果）
if refinement_levels >= 1
     N = N_list(end);
     h = (2*pi)/N;
     vec_p = linspace(0, 2*pi, N+1); vec_p = vec_p(1:end-1);
     [X, Y] = meshgrid(vec_p);
     V_ex = v_exact_func(T_final, X, Y);
     W_ex = w_exact_func(T_final, X, Y);
     
     figure;
     subplot(1,2,1); contourf(X, Y, V, 20); title('Numerical v'); axis equal; colorbar;
     subplot(1,2,2); contourf(X, Y, V_ex, 20); title('Exact v'); axis equal; colorbar;
 end