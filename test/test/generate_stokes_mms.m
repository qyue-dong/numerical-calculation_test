%% generate_unsteady_mms.m (Higher Frequency)
% 修正：将频率从 pi 提升到 2*pi，增大空间误差占比
clear; clc;
syms x y t nu alpha real

% 1. 定义非定常精确解 (频率翻倍)
% 使用 sin(2*pi*x) 依然满足 u0 = -u2 的幽灵点假设(在x=0,1处导数为常数，局部近似奇函数)
% 且 sin(2*pi*0)=0, sin(2*pi*1)=0 满足 Dirichlet
u_exact = exp(-t) * sin(2*pi*x) * sin(2*pi*y);
v_exact = exp(-t) * sin(2*pi*x) * sin(2*pi*y); 
p_exact = sin(2*pi*x) * cos(2*pi*y) * exp(-t);

% 2. 反推 u_raw
u_bar = u_exact; v_bar = v_exact;
u_raw = u_bar - alpha^2 * (diff(u_bar,x,2) + diff(u_bar,y,2));
v_raw = v_bar - alpha^2 * (diff(v_bar,x,2) + diff(v_bar,y,2));

% 3. 计算各项
dt_u = diff(u_raw, t);
dt_v = diff(v_raw, t);
lap_u = diff(u_raw,x,2) + diff(u_raw,y,2);
lap_v = diff(v_raw,x,2) + diff(v_raw,y,2);

% 解析对流项
adv_u = u_bar * diff(u_raw, x) + v_bar * diff(u_raw, y);
adv_v = u_bar * diff(v_raw, x) + v_bar * diff(v_raw, y);

dp_dx = diff(p_exact, x);
dp_dy = diff(p_exact, y);

% 总源项
f_x = dt_u + adv_u - nu * lap_u + dp_dx;
f_y = dt_v + adv_v - nu * lap_v + dp_dy;

fprintf('正在生成高频源项 mms_source_func.m ...\n');
matlabFunction(f_x, f_y, dp_dx, dp_dy, adv_u, adv_v, u_raw, v_raw, u_bar, v_bar, ...
    'File', 'mms_source_func', 'Vars', {t, x, y, nu, alpha});
fprintf('生成完毕！\n');