% 版本: v1.0 - 经过验证的线性斯托克斯求解器

function [V_next, W_next] = Your_NS_Alpha_Solver(V_prev, W_prev, ...
                                               h, tau, nu, Sv_matrix, Sw_matrix, ops)

% 这一行 'function ...' 必须是文件的第一行，前面不能有任何东西。
    
    N = size(V_prev, 1);
    N2 = N^2;

    % --- 向量化 ---
    v_prev = reshape(V_prev, N2, 1);
    w_prev = reshape(W_prev, N2, 1);
    Sv = reshape(Sv_matrix, N2, 1);
    Sw = reshape(Sw_matrix, N2, 1);

    % --- a. 预测阶段 ---
    LHS = speye(N2) - 0.5 * tau * nu * ops.L_vel;
    RHS_v = (speye(N2) + 0.5 * tau * nu * ops.L_vel) * v_prev + tau * Sv;
    RHS_w = (speye(N2) + 0.5 * tau * nu * ops.L_vel) * w_prev + tau * Sw;
    
    v_star = LHS \ RHS_v;
    w_star = LHS \ RHS_w;

    % --- b. 压力泊松方程 ---
    div_u_star = ops.Div_u * v_star + ops.Div_v * w_star;
    RHS_p = (1/tau) * div_u_star;
    RHS_p = RHS_p - mean(RHS_p);
    
    p_next = ops.Lap_p \ RHS_p;
    p_next = p_next - mean(p_next);

    % --- c. 速度校正 ---
    v_next_vec = v_star - tau * (ops.Grad_x * p_next);
    w_next_vec = w_star - tau * (ops.Grad_y * p_next);

    % --- d. 反向量化 ---
    V_next = reshape(v_next_vec, N, N);
    W_next = reshape(w_next_vec, N, N);

end % 这个 end 对应于 function```