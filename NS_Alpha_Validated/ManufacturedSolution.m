classdef ManufacturedSolution
    properties
        nu, alpha, U0
    end
    methods
        function obj = ManufacturedSolution(nu_val, alpha_val)
            obj.nu = nu_val; obj.alpha = alpha_val; obj.U0 = 1.0;
        end
        function [v, w, p, v_bar, w_bar] = getSolution(obj, t, x, y)
            decay_vel = exp(-2 * obj.nu * t);
            decay_p = exp(-4 * obj.nu * t);
            v = -obj.U0 * cos(x) .* sin(y) * decay_vel;
            w =  obj.U0 * sin(x) .* cos(y) * decay_vel;
            p = - (obj.U0^2 / 4) * (cos(2*x) + cos(2*y)) * decay_p;
            lambda = 2; helmholtz_factor = 1 / (1 + obj.alpha^2 * lambda);
            v_bar = v * helmholtz_factor; w_bar = w * helmholtz_factor;
        end
      
% --- 最终正确的代码 ---
function [Sv, Sw] = getSourceTerms(obj, t, x, y)
    % 这个函数现在计算的是 *线性斯托克斯方程* 的源项。
    
    [v, w, p, ~, ~] = obj.getSolution(t, x, y); % 不再需要 v_bar, w_bar
    decay_p = exp(-4 * obj.nu * t);

    % 时间导数
    v_t = -2 * obj.nu * v;
    w_t = -2 * obj.nu * w;
    
    % 拉普拉斯项
    delta_v = -2 * v;
    delta_w = -2 * w;
    
    % 压力梯度
    p_x = - (obj.U0^2 / 4) * (-2 * sin(2*x)) * decay_p;
    p_y = - (obj.U0^2 / 4) * (-2 * sin(2*y)) * decay_p;

    % 源项 S = ∂u/∂t - νΔu + ∇p
    Sv = v_t - obj.nu * delta_v + p_x; % <--- 移除了非线性项
    Sw = w_t - obj.nu * delta_w + p_y; % <--- 移除了非线性项
end


        function [v0, w0] = getInitialCondition(obj, x, y)
            [v0, w0, ~, ~, ~] = obj.getSolution(0, x, y);
        end
    end
end