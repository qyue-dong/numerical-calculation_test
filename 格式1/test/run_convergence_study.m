function run_convergence_study()
%%==============================================================
%  NS-ALPHA MMS CONVERGENCE TEST (FULLY DISCRETE-CONSISTENT)
%===============================================================

clear; clc; close all;

Re = 100;
nu = 1/Re;
L = 1;
T_final = 0.1;

grid_sizes = [16 32 64 128];
errors = zeros(length(grid_sizes),1);

fprintf("=============================================================\n");
fprintf("  NS-alpha (Compact Scheme) MMS Convergence Test\n");
fprintf("=============================================================\n");

for k = 1:length(grid_sizes)

    N = grid_sizes(k);
    h = L/(N-1);
    dt = 0.5*h^2;
    alpha = h;

    fprintf("Grid N=%d, h=%.4f, dt=%.2e\n",N,h,dt);

    % ----------------------------------------------------------
    % Grid + operators
    % ----------------------------------------------------------
    x = linspace(0,1,N);
    [X,Y] = meshgrid(x,x);

    ops = build_spatial_operators(N,h);

    I = speye(N^2);

    H1 = ops.H1;
    H2 = ops.H2;
    Dx = ops.A1;
    Dy = ops.A2;

    Lp = ops.Lp_N;
    Lp(1,:)=0; Lp(:,1)=0; Lp(1,1)=1;

    dec_p = decomposition(Lp,'lu');

    % Velocity matrices
    M_L = H1/dt - 0.5*nu*H2;
    M_R = H1/dt + 0.5*nu*H2;

    % Boundary mask
    bd = false(N,N);
    bd(1,:)=1; bd(end,:)=1; bd(:,1)=1; bd(:,end)=1;
    idx_bd = find(bd);

    M_L = apply_bc(M_L, idx_bd);

    dec_vel = decomposition(M_L,'lu');

    % ----------------------------------------------------------
    % Initial condition
    % ----------------------------------------------------------
    t = 0;
    [~,~,u,v,ub,vb,~] = mms_source_func_discrete(t,X(:),Y(:),nu,alpha,ops);

    % ----------------------------------------------------------
    % Time stepping
    % ----------------------------------------------------------
    nsteps = ceil(T_final/dt);

    for n = 1:nsteps

        t_next = t + dt;

        % MMS forcing at t+dt/2
        [fx,fy,~,~,~,~,p_ex] = ...
            mms_source_func_discrete(t+0.5*dt,X(:),Y(:),nu,alpha,ops);

        Fx = H1*fx;
        Fy = H1*fy;

        % Advection
        adv_u = -( ub .* (Dx*u) + vb .* (Dy*u) );
        adv_v = -( ub .* (Dx*v) + vb .* (Dy*v) );

        % RHS
        rhs_u = M_R*u + dt*(adv_u + Fx);
        rhs_v = M_R*v + dt*(adv_v + Fy);

        % Exact BC at t_next
        [~,~,u_ex_b,v_ex_b,~,~,~] = ...
            mms_source_func_discrete(t_next,X(:),Y(:),nu,alpha,ops);

        rhs_u(idx_bd) = u_ex_b(idx_bd);
        rhs_v(idx_bd) = v_ex_b(idx_bd);

        % Predictor
        u_star = dec_vel \ rhs_u;
        v_star = dec_vel \ rhs_v;

        % Pressure correction
        div_star = Dx*u_star + Dy*v_star;
        rhs_p = div_star/dt;
        rhs_p(1) = 0;

        p_new = dec_p \ rhs_p;

        % Correct velocity
        u_new = u_star - dt*(ops.Grad_x*p_new);
        v_new = v_star - dt*(ops.Grad_y*p_new);

        % enforce BC
        u_new(idx_bd) = u_ex_b(idx_bd);
        v_new(idx_bd) = v_ex_b(idx_bd);

        % Filtered velocity = exact in MMS
        ub = u_new;
        vb = v_new;

        u = u_new;
        v = v_new;
        t = t_next;

    end

    % ----------------------------------------------------------
    % Error
    % ----------------------------------------------------------
    [~,~,u_ex,v_ex,~,~,~] = ...
        mms_source_func_discrete(T_final,X(:),Y(:),nu,alpha,ops);

    err_u = norm(u-u_ex)/norm(u_ex);
    err_v = norm(v-v_ex)/norm(v_ex);

    total = sqrt(err_u^2 + err_v^2);

    errors(k) = total;

    fprintf("  L2 error = %.4e\n", total);

end

% --------------------------------------------------------------
% Convergence
% --------------------------------------------------------------
hvals = L./(grid_sizes-1);
rates = log2(errors(1:end-1)./errors(2:end));

fprintf("\n--- Convergence Rates ---\n");
for i=1:length(rates)
    fprintf("%d -> %d : %.4f\n", grid_sizes(i),grid_sizes(i+1),rates(i));
end

% Plot
figure;
loglog(hvals,errors,'bo-','LineWidth',2,'MarkerSize',8);
hold on
loglog(hvals, errors(1)*(hvals/hvals(1)).^4, 'k--','LineWidth',1.5);
loglog(hvals, errors(1)*(hvals/hvals(1)).^3, 'r:','LineWidth',1.5);
xlabel('Grid spacing h');
ylabel('L2 error');
legend('Error','4th order','3rd order','Location','best');
grid on;

end

function A = apply_bc(A,idx)
for k = idx(:)'
    A(k,:) = 0;
    A(k,k) = 1;
end
end
