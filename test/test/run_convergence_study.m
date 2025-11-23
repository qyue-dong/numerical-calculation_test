function run_convergence_study()
%% Final Verification: Euler + h^4 Scaling (Bulletproof 4th Order)
% 策略：
%   1. 求解器：全隐式欧拉 (Implicit Euler) -> 消除 CN 震荡，保证绝对稳定
%   2. 时间步长：dt = 5.0 * h^4 -> 强制时间误差与空间误差同阶 (4阶)
%   3. 模拟时长：极短 -> 保证计算在几秒内完成

clear; clc; close all;

%% 1. 实验设置
Re = 100;
nu = 1/Re;
L = 1.0;
% 只需要跑很短的时间就能看清空间精度
% 因为 h^4 极小，跑太久没有意义且太慢
T_final = 1e-3; 

grid_sizes = [12, 16, 20, 24, 32]; 
errors_u = zeros(length(grid_sizes), 1);

fprintf('============================================================\n');
fprintf('Final Verification: Implicit Euler with dt ~ h^4\n');
fprintf('============================================================\n');

for k = 1:length(grid_sizes)
    N = grid_sizes(k);
    h = L / (N-1); 
    alpha = h;
    
    % --- 核心策略：dt ~ h^4 ---
    % 隐式欧拉是一阶时间精度 Error ~ O(dt)
    % 我们要验证四阶空间精度 Error ~ O(h^4)
    % 令 dt = C * h^4，则总误差 Error ~ O(h^4)
    dt = 2.0 * h^4; 
    
    % --- A. 构建算子 ---
    ops = build_spatial_operators(N, h);
    x_1d = linspace(0, L, N);
    [X, Y] = meshgrid(x_1d, x_1d);
    
    % --- B. 提取算子 ---
    I = speye(N^2);
    % 4阶拉普拉斯 (Ghost Point 修正版)
    Lp_4th = ops.A3 + ops.A5; 
    
    % --- C. 组装时间步进矩阵 (Implicit Euler) ---
    % (I - dt*nu*Lp) * u_new = u_old + dt*Source
    M_LHS = I - dt * nu * Lp_4th;
    M_RHS_Mass = I;
    
    % 滤波算子
    M_filter = I - alpha^2 * Lp_4th;
    
    % --- D. 边界处理 (Lifting) ---
    bd_mask = false(N, N);
    bd_mask(1,:) = true; bd_mask(end,:) = true;
    bd_mask(:,1) = true; bd_mask(:,end) = true;
    idx_bd = find(bd_mask);
    
    % 预分解
    [M_LHS_Mod, ~] = apply_bc_lifting(M_LHS, zeros(N^2,1), idx_bd, 0); 
    dA_vel = decomposition(M_LHS_Mod, 'lu');
    
    [M_filter_Mod, ~] = apply_bc_lifting(M_filter, zeros(N^2,1), idx_bd, 0);
    dA_filter = decomposition(M_filter_Mod, 'lu');
    
    % --- E. 初始化 ---
    t = 0;
    [~, ~, ~, ~, ~, ~, u, v, ~, ~] = mms_source_func(t, X(:), Y(:), nu, alpha);
    
    % --- F. 时间循环 ---
    num_steps = ceil(T_final / dt);
    fprintf('Grid %dx%d (dt=%.1e, Steps=%d)... ', N, N, dt, num_steps);
    
    for n = 1:num_steps
        t_next = t + dt;
        
        % 1. 计算源项 (t_next)
        % Implicit Euler 使用 t_next 时刻的源项
        [fx, fy, dpx, dpy, adv_u_ana, adv_v_ana, u_ex, v_ex, ~, ~] ...
            = mms_source_func(t_next, X(:), Y(:), nu, alpha);
        
        % 2. 构造纯净源项 (Source = f - grad_p - adv)
        Source_u = fx - dpx - adv_u_ana;
        Source_v = fy - dpy - adv_v_ana;
        
        % 3. 组装 RHS (Euler: u_old + dt * Source)
        b_u = u + dt * Source_u;
        b_v = v + dt * Source_v;
        
        % 4. 应用边界 (Lifting)
        [~, b_u_mod] = apply_bc_lifting(M_LHS, b_u, idx_bd, u_ex(idx_bd));
        [~, b_v_mod] = apply_bc_lifting(M_LHS, b_v, idx_bd, v_ex(idx_bd));
        
        % 5. 求解
        u_new = dA_vel \ b_u_mod;
        v_new = dA_vel \ b_v_mod;
        
        % 6. 滤波 (后处理，不影响 u 的收敛性)
        rhs_ub = u_new; rhs_vb = v_new;
        [~, rhs_ub] = apply_bc_lifting(M_filter, rhs_ub, idx_bd, 0); % 简化的滤波边界
        [~, rhs_vb] = apply_bc_lifting(M_filter, rhs_vb, idx_bd, 0);
        u_bar = dA_filter \ rhs_ub;
        v_bar = dA_filter \ rhs_vb;
        
        % 更新
        u = u_new; v = v_new;
        t = t_next;
        
        if any(isnan(u)), error('NaN'); end
    end
    
    % --- G. 误差计算 ---
    [~, ~, ~, ~, ~, ~, u_final, v_final, ~, ~] = mms_source_func(t, X(:), Y(:), nu, alpha);
    err = sqrt(sum((u - u_final).^2 + (v - v_final).^2) / length(u));
    errors_u(k) = err;
    
    fprintf('Error = %.4e\n', err);
end

%% 绘图
h_vals = L ./ (grid_sizes - 1);
figure('Color','w');
loglog(h_vals, errors_u, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;

% 绘制参考线
ref_h = h_vals;
% 调整参考线位置，方便对比
start_err = errors_u(1);
loglog(ref_h, start_err*(ref_h/ref_h(1)).^4, 'k--', 'LineWidth', 1.5, 'DisplayName', '4th Order');
loglog(ref_h, start_err*(ref_h/ref_h(1)).^3, 'r:', 'LineWidth', 1.5, 'DisplayName', '3rd Order');

legend('Our Simulation', '4th Order', '3rd Order', 'Location', 'best');
grid on; xlabel('h'); ylabel('L2 Error');
title('Final Proof: Implicit Euler (dt ~ h^4)');

% 通用收敛阶计算公式
rates = log(errors_u(1:end-1)./errors_u(2:end)) ./ log(h_vals(1:end-1)./h_vals(2:end))';
fprintf('\n--- Final Convergence Rates ---\n');
disp(rates);

end

% --- Lifting 边界处理 ---
function [A_mod, b_mod] = apply_bc_lifting(A, b, idx_bd, u_bd_vals)
    [n, ~] = size(A);
    u_known = zeros(n, 1);
    u_known(idx_bd) = u_bd_vals;
    b_mod = b - A * u_known; 
    I = speye(n);
    d = ones(n, 1); d(idx_bd) = 0;
    D = spdiags(d, 0, n, n);
    A_mod = D * A + (I - D); 
    b_mod(idx_bd) = u_bd_vals; 
end