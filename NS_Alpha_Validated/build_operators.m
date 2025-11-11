% --- build_operators.m ---
% 版本: v1.0 - 经过验证的纯二阶稳定版
function ops = build_operators(N, h)
    fprintf('>>> 正在运行【v1.0 - 纯二阶稳定版】build_operators <<<\n');
    J = N;
    I_J = speye(J);
    e = ones(J, 1);
    
    % --- 1D 核心矩阵 ---
    C1 = spdiags([-e, e], [-1, 1], J, J);
    C1(1, J) = -1; C1(J, 1) = 1;
    C2 = spdiags([e, -2*e, e], [-1, 0, 1], J, J);
    C2(1, J) = 1; C2(J, 1) = 1;

    % --- 2D 基础算子 ---
    ops.A1 = (1/(2*h)) * kron(I_J, C1); % d/dx
    ops.A2 = (1/(2*h)) * kron(C1, I_J); % d/dy
    ops.A3 = (1/h^2) * kron(I_J, C2); % d^2/dx^2
    ops.A5 = (1/h^2) * kron(C2, I_J); % d^2/dy^2

    % --- 速度方程专用算子 ---
    ops.L_vel = ops.A3 + ops.A5; % 速度拉普拉斯 (Δ)
    
    % --- 压力-速度耦合专用算子 ---
    ops.Grad_x = ops.A1;
    ops.Grad_y = ops.A2;
    ops.Div_u = ops.A1;
    ops.Div_v = ops.A2;
    ops.Lap_p = ops.Grad_x * ops.Grad_x + ops.Grad_y * ops.Grad_y;
    
    fprintf('所有【纯二阶稳定版】算子构建完成。\n');
end