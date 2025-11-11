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
        function [Sv, Sw] = getSourceTerms(obj, t, x, y)
            [v, w, ~, v_bar, w_bar] = obj.getSolution(t, x, y);
            decay_vel = exp(-2 * obj.nu * t); decay_p = exp(-4 * obj.nu * t);
            v_t = -2 * obj.nu * v; w_t = -2 * obj.nu * w;
            delta_v = -2 * v; delta_w = -2 * w;
            p_x = - (obj.U0^2 / 4) * (-2 * sin(2*x)) * decay_p;
            p_y = - (obj.U0^2 / 4) * (-2 * sin(2*y)) * decay_p;
            w_x =  obj.U0 * cos(x) .* cos(y) * decay_vel;
            v_y = -obj.U0 * cos(x) .* cos(y) * decay_vel;
            omega_z = w_x - v_y;
            Sv = v_t - obj.nu * delta_v - w_bar .* omega_z + p_x;
            Sw = w_t - obj.nu * delta_w + v_bar .* omega_z + p_y;
        end
        function [v0, w0] = getInitialCondition(obj, x, y)
            [v0, w0, ~, ~, ~] = obj.getSolution(0, x, y);
        end
    end
end