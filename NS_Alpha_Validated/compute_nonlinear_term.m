function [NL_v, NL_w] = compute_nonlinear_term(V, W, V_bar, W_bar, ops)
    % compute_nonlinear_term: 计算NS-alpha模型中的非线性项。
    % 非线性项 N(u) = ((∇×u)×ū)
    % 其x分量: -w_bar * (∂w/∂x - ∂v/∂y)
    % 其y分量:  v_bar * (∂w/∂x - ∂v/∂y)
    %
    % 输入:
    %   V, W, V_bar, W_bar: 当前时刻的速度场 (N x N 矩阵)
    %   ops: 算子结构体
    % 输出:
    %   NL_v, NL_w: 非线性项的x和y分量 (N^2 x 1 列向量)

    N2 = size(V, 1) * size(V, 2);

    % --- 向量化输入 ---
    v = reshape(V, N2, 1);
    w = reshape(W, N2, 1);
    v_bar = reshape(V_bar, N2, 1);
    w_bar = reshape(W_bar, N2, 1);

    % --- 计算涡度 (ω_z = ∂w/∂x - ∂v/∂y) ---
    % 我们使用与压力系统兼容的二阶算子 Grad_x, Grad_y
   
   w_x = ops.Grad_x * w; % <--- 正确！
   v_y = ops.Grad_y * v; % <--- 正确！
   omega_z = w_x - v_y;
   
    % --- 计算非线性项 ---
    % 这里的乘法是逐点相乘 (element-wise)
    NL_v = -w_bar .* omega_z;
    NL_w =  v_bar .* omega_z;
end