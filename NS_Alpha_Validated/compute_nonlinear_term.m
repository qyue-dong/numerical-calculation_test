% 版本: v1.1 - 最终统一命名版
function [NL_v, NL_w] = compute_nonlinear_term(V, W, V_bar, W_bar, ops)
    N2 = size(V, 1) * size(V, 2);
    v = reshape(V, N2, 1);
    w = reshape(W, N2, 1);
    v_bar = reshape(V_bar, N2, 1);
    w_bar = reshape(W_bar, N2, 1);

    % --- 计算涡度 (ω_z = ∂w/∂x - ∂v/∂y) ---
    w_x = ops.Dx * w; % << 使用新的、清晰的算子名
    v_y = ops.Dy * v; % << 使用新的、清晰的算子名
    omega_z = w_x - v_y;

    % --- 计算非线性项 ---
    NL_v = -w_bar .* omega_z;
    NL_w =  v_bar .* omega_z;
end