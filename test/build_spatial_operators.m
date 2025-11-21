function ops = build_spatial_operators(N,h)

% ===== circulant tri-diagonal constructor =====
circ = @(a,b,c,N) spdiags([c*ones(N,1), b*ones(N,1), a*ones(N,1)],[-1 0 1],N,N) ...
                   + sparse(1,N,a) + sparse(N,1,c);

I = speye(N);

% ===============================
% C1 (first derivative circulant)
% ===============================
C1 = circ(+1,0,-1,N);   % matches your screenshot for C1
A1 = kron(C1, I) / (2*h);
A2 = kron(I, C1) / (2*h);

% ===============================
% C2 (Laplacian core circulant)
% ===============================
C2 = circ(1,-2,1,N);
A3 = kron(C2, I) / (h^2);
A5 = kron(I, C2) / (h^2);

% ===============================
% A4 = (h^2/12)*A5  (from screenshot)
% ===============================
A4 = (h^2/12)*A5;

% ===============================
% C3, C4, C5 (higher-order stencils)
% ===============================

% your screenshot says:
% C3: diag = 10, offdiag = 1
C3 = circ(1,10,1,N);

% C4: diag = 14, offdiag = -1
C4 = circ(-1,14,-1,N);

% C5: diag = 8, offdiag = -1
C5 = circ(-1,8,-1,N);

% build B operators exactly as doc:
B1 = kron(C3,I) / 12;
B2 = kron(C4,I) / 12;
B3 = kron(I,C4) / 12;
B4 = kron(C5,I) / 6;
B5 = kron(I,C5) / 6;

% ===============================
% H1, H2 definitions
% ===============================
H1 = B1 + A4;
H2 = B2*A3 + B3*A5;

% ===============================
% pressure gradient / Laplacian
% ===============================
Grad_x = B4*A1;
Grad_y = B5*A2;

% pressure Laplace (Neumann)
Lp_N = A3 + A5;

ops.A1 = A1; ops.A2 = A2;
ops.A3 = A3; ops.A5 = A5;
ops.A4 = A4;

ops.B1 = B1; ops.B2 = B2; ops.B3 = B3;
ops.B4 = B4; ops.B5 = B5;

ops.H1 = H1; ops.H2 = H2;

ops.Grad_x = Grad_x;
ops.Grad_y = Grad_y;

ops.Lp_N   = Lp_N;

end
