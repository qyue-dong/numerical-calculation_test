function [fx, fy, u_ex, v_ex, ub_ex, vb_ex, p_ex] = ...
    mms_source_func_discrete(t, X, Y, nu, alpha, ops)
% =====================================================================
%  MMS FOR THE NS-ALPHA COMPACT SCHEME (DISCRETE CONSISTENT)
%  Fully discrete:
%     H1*u_t + H1*(u_bar·∇u) - ν H2*u + ∇p = f
% =====================================================================

% 1. Exact solution -----------------------------------------------------
u_ex = exp(-nu*t) .* (sin(pi*X).^2 .* sin(pi*Y) .* cos(pi*Y));
v_ex = exp(-nu*t) .* (sin(pi*Y).^2 .* sin(pi*X) .* cos(pi*X));
p_ex = exp(-nu*t) .* (cos(pi*X) .* sin(pi*Y));

ub_ex = u_ex;      % MMS chooses u_bar = u
vb_ex = v_ex;

% Flatten
u  = u_ex(:);
v  = v_ex(:);
ub = ub_ex(:);
vb = vb_ex(:);
p  = p_ex(:);

% 2. Exact time derivative ----------------------------------------------
u_t = -nu * u;
v_t = -nu * v;

% 3. Discrete advection -------------------------------------------------
Du_u = ops.A1 * u;
Dv_u = ops.A2 * u;
Du_v = ops.A1 * v;
Dv_v = ops.A2 * v;

adv_u = ub .* Du_u + vb .* Dv_u;
adv_v = ub .* Du_v + vb .* Dv_v;

% 4. Discrete Laplacian (H2) --------------------------------------------
L_u = ops.H2 * u;
L_v = ops.H2 * v;

% 5. Discrete pressure gradient -----------------------------------------
dpdx = ops.Grad_x * p;
dpdy = ops.Grad_y * p;

% 6. Final MMS forcing (fully discrete) ---------------------------------
Fx = ops.H1 * (u_t + adv_u) - nu*L_u + dpdx;
Fy = ops.H1 * (v_t + adv_v) - nu*L_v + dpdy;

fx = Fx;
fy = Fy;

end
