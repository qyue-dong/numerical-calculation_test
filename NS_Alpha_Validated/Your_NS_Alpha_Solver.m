% --- Your_NS_Alpha_Solver.m ---
% 版本: v1.1 - 最终统一命名版
function [V_next, W_next] = Your_NS_Alpha_Solver(V_prev, W_prev, ...
                                               h, tau, nu, Sv_matrix, Sw_matrix, ops)
    N = size(V_prev, 1);
    N2 = N^2;

    v_prev = reshape(V_prev, N2, 1);
    w_prev = reshape(W_prev, N2, 1);
    Sv = reshape(Sv_matrix, N2, 1);
    Sw = reshape(Sw_matrix, N2, 1);

    % --- a. 预测阶段 ---
    LHS = speye(N2) - 0.5 * tau * nu * ops.Laplacian; % << 使用统一的拉普拉斯算子
    RHS_v = (speye(N2) + 0.5 * tau * nu * ops.Laplacian) * v_prev + tau * Sv;
    RHS_w = (speye(N2) + 0.5 * tau * nu * ops.Laplacian) * w_prev + tau * Sw;
    
    v_star = LHS \ RHS_v;
    w_star = LHS \ RHS_w;

    % --- b. 压力泊松方程 ---
    div_u_star = ops.Dx * v_star + ops.Dy * w_star; % << 使用统一的 Dx, Dy
    RHS_p = (1/tau) * div_u_star;
    RHS_p = RHS_p - mean(RHS_p);
    
    p_next = ops.Laplacian \ RHS_p; % << 使用统一的拉普拉斯算子
    p_next = p_next - mean(p_next);

    % --- c. 速度校正 ---
    v_next_vec = v_star - tau * (ops.Dx * p_next); % << 使用统一的 Dx
    w_next_vec = w_star - tau * (ops.Dy * p_next); % << 使用统一的 Dy

    % --- d. 反向量化 ---
    V_next = reshape(v_next_vec, N, N);
    W_next = reshape(w_next_vec, N, N);
end