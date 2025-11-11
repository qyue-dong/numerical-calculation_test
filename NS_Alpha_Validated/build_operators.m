% 版本: v1.1 - 最终统一命名版
function ops = build_operators(N, h)
    fprintf('>>> 正在运行【v1.1 - 最终统一命名版】build_operators <<<\n');
    J = N;
    I_J = speye(J);
    e = ones(J, 1);
    
    C1 = spdiags([-e, e], [-1, 1], J, J);
    C1(1, J) = -1; C1(J, 1) = 1;
    C2 = spdiags([e, -2*e, e], [-1, 0, 1], J, J);
    C2(1, J) = 1; C2(J, 1) = 1;

    % --- 2D 基础算子 ---
    ops.Dx = (1/(2*h)) * kron(I_J, C1); % d/dx (一阶)
    ops.Dy = (1/(2*h)) * kron(C1, I_J); % d/dy (一阶)
    ops.Dxx = (1/h^2) * kron(I_J, C2);  % d^2/dx^2 (二阶)
    ops.Dyy = (1/h^2) * kron(C2, I_J);  % d^2/dy^2 (二阶)

    % --- 组合算子 ---
    ops.Laplacian = ops.Dxx + ops.Dyy; % 拉普拉斯算子
end