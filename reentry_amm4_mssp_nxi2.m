function reentry_amm4_mssp()
% reentry_amm4_mssp  按 AMM_4 + MSSP_1 谱(PC)混合不确定性方法求解再入滑翔鲁棒轨迹。
%
% 方法来源：
%   MSSP_1（Legendre PC）：mu = y0，sigma^2 = sum_{j>=1} y_j^2 <phi_j^2>；
%   AMM_4（随机 xi + 区间 eta 混合）：
%     区间均值  [mu_F] = beta_0(eta)   （Eq.16）
%     区间方差  [sigma^2_F] = sum_{j>=1} beta_j(eta)^2 <phi_j^2>  （Eq.17）
%     区间上界  sigma^2 = sum_{j>=1} max(beta_j^2, betaBar_j^2) <phi_j^2> （Eq.30）
%   优化模型（AMM_4 Eq.4）：min w1*muBar_f + w2*sigmaBar_f，s.t. muBar_g + w3*sigmaBar_g <= 0。
%
% 实现要点：
%   - xi（均匀随机，n_xi=2：theta0、phi0）用归一化 Legendre 基做 Galerkin 投影，
%     得到谱系数 beta_j(eta)；均值/方差用“系数平方 x <phi_j^2>”的谱形式。
%   - eta（区间，n_eta=3：rho、CL、CD）用 Chebyshev-Lobatto 节点，对 eta 取上下界
%     （区间均值上下界 + 区间方差上界），即 AMM_4 的区间算术/包含函数思想。
%   - 终端：区间均值落在目标邻域（软），区间方差上界（epigraph）进入目标惩罚。
%   - 路径：对每个 eta_q 施加 mu_r >= pc + kC*sigma_r（线性化背向），最坏覆盖 eta。
%
% 需要 YALMIP + MOSEK。运行前把 YALMIP 加入 MATLAB 路径（MOSEK 路径自动添加）。

% 自动添加 MOSEK 路径（若存在）
if isempty(which('mosekopt'))
    msk = 'C:\Program Files\Mosek\11.0\toolbox\r2019b';
    if exist(msk, 'dir'), addpath(msk); end
end

if isempty(which('sdpvar'))
    error('未找到 YALMIP，请先安装并把 YALMIP 目录加入 MATLAB 路径。');
end

%% ===================== 物理参数 =====================
R0    = 6378e3;  g0 = 9.81;  m = 37362;  Sref = 149.4;
omega = 7.2921e-5*sqrt(R0/g0);
kQ = 9.4369e-5;  Qmax = 800e3;  qmax = 14364;  gmax = 2.5;
vmin = 0;  vmax = 165*pi/180;

H0 = 121.9e3;  Hf = 30.48e3;
r0 = (H0+R0)/R0;  rf = (Hf+R0)/R0;
theta0 = 237*pi/180;   thetaf = 279*pi/180;
fai0   = -25*pi/180;   faif   = 28.61*pi/180;
V0     = 7626/sqrt(g0*R0);  Vf = 908.15/sqrt(g0*R0);
gamma0   = -1.2493*pi/180;  gammaf = [-6,0]*pi/180;
polyfai0 = 45*pi/180;       polyfaif = [85,90]*pi/180;

delta   = [20000/R0, 20*pi/180, 20*pi/180, 20*pi/180, 70*pi/180];
epsilon = [200/R0, 0.05*pi/180, 0.05*pi/180, 0.05*pi/180, 0.05*pi/180];

kC = 1;  kJ = 3;     % 终端散布权重（原 3 太小，相对飞行时间可忽略）
N_points = 101;   % 301=论文精度，101=快速
p = 2;
solver_used = 'mosek';   % 首选 MOSEK；许可证失效时 scp_iteration 内自动回退 coneprog

e0 = 1/r0 - V0^2/2;
ef = 1/rf - Vf^2/2;
de = (ef-e0)/(N_points-1);

%% ===================== 混合不确定性配置 =====================
% xi（均匀随机，Legendre）：初始状态不确定量，此处取 theta0、phi0 两个（n_xi=2）。
%   可扩展为全部 5 个初始状态，但全张量配点数会按 N_xi^n_xi 指数增长。
n_xi = 2;   N_xi = 2;                       % Legendre 维数与每维 Gauss-Legendre 点数
sigma_x0 = [0.3*pi/180, 0.3*pi/180];         % 初始状态半宽（此部分开环控制难以消除，取小）

% eta（区间，Chebyshev）：rho、CL、CD（m=3），Chebyshev-Gauss-Lobatto 配点（含端点 ±1）。
n_eta = 3;  N_eta = 3;                       % Chebyshev 维数与每维 Chebyshev-Gauss-Lobatto 点数
delta_eta = [0.10, 0.10, 0.10];              % [rho, CL, CD] 区间半宽（此部分开环控制可抑制，取大）

% ---- 一维求积 ----
[xL, wL] = gauss_legendre(N_xi);
[xC, wC] = cheb_lobatto(N_eta);
[xi_nodes,  w_xi ] = tensor_quad(repmat({xL},1,n_xi), repmat({wL},1,n_xi));
[eta_nodes, w_eta] = tensor_quad(repmat({xC},1,n_eta), repmat({wC},1,n_eta));

Q_xi  = size(xi_nodes,1);
Q_eta = size(eta_nodes,1);
Q     = Q_xi * Q_eta;

% 初始状态随 xi 的变化（只扰动 theta0、phi0，即第 2、3 个状态分量）
X0mat = [r0*ones(Q_xi,1), theta0 + sigma_x0(1)*xi_nodes(:,1), ...
         fai0 + sigma_x0(2)*xi_nodes(:,2), gamma0*ones(Q_xi,1), polyfai0*ones(Q_xi,1)];

% eta -> 动力学参数扰动 [Delta rho, Delta CL, Delta CD]
Seta = [delta_eta(1)*eta_nodes(:,1), delta_eta(2)*eta_nodes(:,2), delta_eta(3)*eta_nodes(:,3)];

% 全配点（q 为主序：先按 eta_q 分组，组内按 xi_s）
SfullR = zeros(Q,3);  X0fullR = zeros(Q,5);  w_fullR = zeros(Q,1);
for q = 1:Q_eta
    idx = (q-1)*Q_xi + (1:Q_xi);
    SfullR(idx,:)  = repmat(Seta(q,:), Q_xi, 1);
    X0fullR(idx,:) = X0mat;
    w_fullR(idx)   = w_eta(q) * w_xi;
end

% ---- 谱(PCE)系数投影：归一化 Legendre 张量基 phi_{k1,k2}=P_k1(xi1)P_k2(xi2) ----
% S=Q_xi（N_xi=2,n_xi=2 => 2x2=4 个完全基函数，Galerkin 投影精确）
Phi_xi = zeros(Q_xi, Q_xi);
ii = 1;
for k1 = 0:1
    for k2 = 0:1
        Phi_xi(:, ii) = legpoly(k1, xi_nodes(:,1)) .* legpoly(k2, xi_nodes(:,2));
        ii = ii + 1;
    end
end
% Galerkin 系数映射（<phi_j^2>=1）：beta_j = C_xiR * f_q，f_q 为 eta_q 处对 xi 的节点值
C_xiR = (Phi_xi .* w_xi.').';         % S x Q_xi
% 谱均值映射（beta_0 行）：对每个 eta_q 对 xi 取期望
A_muR = kron(eye(Q_eta), C_xiR(1,:));  % Q_eta x Q，C_xiR(1,:)=w_xi（phi_0=1）
% 谱方差二次型（MSSP_1 Eq.10 / AMM_4 Eq.17 的节点矩阵形式）：
%   Var[f]=sum_{j>=1} beta_j^2 = f'*(C_xiR(2:end,:)'*C_xiR(2:end,:))*f
M_xiR = C_xiR(2:end,:)' * C_xiR(2:end,:);   % 恒等于 diag(w_xi)-w_xi*w_xi'（Parseval）

% ---- 标称 eta 块（eta=0，Chebyshev-Lobatto N_eta=3 含 0 节点），性能项在此对 xi 求期望 ----
[~, q0] = min(sum(Seta.^2, 2));      % 最接近 eta=0 的节点下标
idx0 = (q0-1)*Q_xi + (1:Q_xi);       % 标称 eta 对应的全配点索引
w_perf_R = zeros(Q, 1);  w_perf_R(idx0) = w_xi;  % 性能权重：仅在标称 eta 处对 xi 求期望

%% ===================== 公共共享状态（嵌套函数使用） =====================
e_values = []; dei = []; dei_1 = [];
% 确定性（标称）阶段：Q=1
Scur   = zeros(1,3);
X0full = [r0, theta0, fai0, gamma0, polyfai0];
w_full = 1;  A_mu = 1;  M_xi = 0;  w_perf = 1;
Qcur = 1;  Qxi = 1;  Qeta = 1;  tr_scale = 1;

init_e_values();
pc = path_constrain_init(zeros(1,3), 1);

sigma_nom = sigma_e_function();
[rk, thetak, faik, gammak, polyfaik] = init_trajectories(sigma_nom);
u1 = cosd(sigma_nom(e_values));  u2 = sind(sigma_nom(e_values));

fprintf('===== 确定性轨迹优化（DO）=====\n');
tic
prev_maxdX = Inf;
for iter = 1:100
    [ok, r_n, t_n, f_n, g_n, p_n, u1n, u2n] = scp_iteration(rk, thetak, faik, gammak, polyfaik, pc, u1, u2);
    if ~ok
        tr_scale = max(tr_scale*0.5, 0.02);
        if tr_scale <= 0.021, fprintf('确定性阶段第 %d 步失败且信赖域耗尽\n', iter); break; end
        continue;
    end
    dX = max_change({rk,thetak,faik,gammak,polyfaik}, {r_n,t_n,f_n,g_n,p_n});
    maxdX = max(dX);
    fprintf('DO iter %2d: max dX = %.3e\n', iter, maxdX);
    rk=r_n; thetak=t_n; faik=f_n; gammak=g_n; polyfaik=p_n; u1=u1n; u2=u2n;
    if maxdX > prev_maxdX, tr_scale = max(tr_scale*0.5, 0.02); else, tr_scale = min(tr_scale*1.15, 1.0); end
    prev_maxdX = maxdX;
    if all(dX < epsilon), break; end
end
DO_X = {rk, thetak, faik, gammak, polyfaik};  DO_u = {u1, u2};
toc
%% ===================== 鲁棒（混合不确定）阶段 =====================
tic
Scur = SfullR;  X0full = X0fullR;  w_full = w_fullR;  A_mu = A_muR;  M_xi = M_xiR;  w_perf = w_perf_R;
Qcur = Q;  Qxi = Q_xi;  Qeta = Q_eta;
pc = path_constrain_init(Seta, Q_eta);

sigma_rob = @(e) atan2d(interp1(e_values, DO_u{2}, e, 'linear','extrap'), ...
                        interp1(e_values, DO_u{1}, e, 'linear','extrap'));
[rk, thetak, faik, gammak, polyfaik] = init_trajectories(sigma_rob);
u1 = DO_u{1};  u2 = DO_u{2};

fprintf('===== 鲁棒轨迹优化（混合 Legendre+Chebyshev + SCO）=====\n');
prev_maxdX = Inf;
for iter = 1:100
    [ok, r_n, t_n, f_n, g_n, p_n, u1n, u2n] = scp_iteration(rk, thetak, faik, gammak, polyfaik, pc, u1, u2);
    if ~ok
        tr_scale = max(tr_scale*0.5, 0.02);
        if tr_scale <= 0.021, fprintf('鲁棒阶段第 %d 步失败且信赖域耗尽\n', iter); break; end
        continue;
    end
    dX = max_change({rk,thetak,faik,gammak,polyfaik}, {r_n,t_n,f_n,g_n,p_n});
    maxdX = max(dX);
    fprintf('ROB iter %2d: max dX = %.3e  (tr=%.3g)\n', iter, maxdX, tr_scale);
    rk=r_n; thetak=t_n; faik=f_n; gammak=g_n; polyfaik=p_n; u1=u1n; u2=u2n;
    if maxdX > prev_maxdX, tr_scale = max(tr_scale*0.5, 0.02); else, tr_scale = min(tr_scale*1.15, 1.0); end
    prev_maxdX = maxdX;
    if all(dX < epsilon), break; end
end
ROB_X = {rk, thetak, faik, gammak, polyfaik};  ROB_u = {u1, u2};
toc
%% ===================== 绘图与验证 =====================
q_nom = find(all(abs(Seta)<1e-12, 2), 1);
if isempty(q_nom), [~, q_nom] = min(sum(Seta.^2,2)); end
idx_nom = (q_nom-1)*Qxi + (1:Qxi);
r_rob = ROB_X{1}(:, idx_nom) * w_xi;   % eta=0 处对 xi 求期望的名义鲁棒轨迹
pc_nom = pc(:, q_nom);

figure('Name','Velocity-Altitude','Color','w');
plot(V_(e_values,DO_X{1}(:,1))*sqrt(g0*R0)/1000, (DO_X{1}(:,1)-1)*R0/1000,'b-','LineWidth',2); hold on;
plot(V_(e_values,r_rob)*sqrt(g0*R0)/1000, (r_rob-1)*R0/1000,'r--','LineWidth',2);
plot(V_(e_values,pc_nom)*sqrt(g0*R0)/1000, (pc_nom-1)*R0/1000,'k-.','LineWidth',1.5);
xlabel('Velocity (km/s)'); ylabel('Altitude (km)'); legend('DO','Hybrid robust','Path constraint','Location','best');
set(gca,'XDir','reverse'); grid on; drawnow;

figure('Name','Bank Angle','Color','w');
plot(e_values, atan2d(DO_u{2},DO_u{1}),'b-','LineWidth',2); hold on;
plot(e_values, atan2d(ROB_u{2},ROB_u{1}),'r--','LineWidth',2);
xlabel('Energy e (non-dim)'); ylabel('Bank angle (deg)'); legend('DO','Hybrid robust','Location','best'); grid on; drawnow;

figure('Name','Control Magnitude','Color','w');
plot(e_values, sqrt(DO_u{1}.^2 + DO_u{2}.^2),'b-','LineWidth',2); hold on;
plot(e_values, sqrt(ROB_u{1}.^2 + ROB_u{2}.^2),'r--','LineWidth',2);
yline(1,'k-.','LineWidth',1.5);
xlabel('Energy e (non-dim)'); ylabel('Control magnitude |u|'); legend('DO','Hybrid robust','|u|=1','Location','best'); grid on; drawnow;

% 所有配点（全张量）轨迹：状态量、控制量、过程约束
plot_all_trajectories(ROB_X, ROB_u);

rng('default');  N_mc = 300;
[t_do,f_do,v_do,MC_DO] = mc_verify(DO_u{1},DO_u{2},N_mc);
[t_ro,f_ro,v_ro,MC_ROB] = mc_verify(ROB_u{1},ROB_u{2},N_mc);
fprintf('\n===== 蒙特卡洛验证（N=%d）=====\n', N_mc);
fprintf('DO  : theta std=%.4f deg, phi std=%.4f deg, 约束违反率=%.2f%%\n', std(t_do)*180/pi, std(f_do)*180/pi, 100*mean(v_do));
fprintf('ROB : theta std=%.4f deg, phi std=%.4f deg, 约束违反率=%.2f%%\n', std(t_ro)*180/pi, std(f_ro)*180/pi, 100*mean(v_ro));

% DO 与 ROB 的标称轨迹、三维打靶轨迹、状态曲线、路径约束、终端散布对比
plot_mc_comparison(DO_X, DO_u, ROB_X, ROB_u, MC_DO, MC_ROB);
fprintf('\n求解完成。\n');

% 导出关键数据到基础工作区，便于运行后直接查看/分析
assignin('base','e_values',e_values);
assignin('base','DO_X',DO_X);        assignin('base','DO_u',DO_u);
assignin('base','ROB_X',ROB_X);      assignin('base','ROB_u',ROB_u);
assignin('base','MC_DO',MC_DO);      assignin('base','MC_ROB',MC_ROB);
assignin('base','MC_DO_terminal',[t_do,f_do,v_do]);
assignin('base','MC_ROB_terminal',[t_ro,f_ro,v_ro]);
assignin('base','Scur',Scur);        assignin('base','X0full',X0full);
assignin('base','Seta',Seta);        assignin('base','w_full',w_full);
assignin('base','Qxi',Qxi);          assignin('base','Qeta',Qeta);

%% ===================== 嵌套函数 =====================

function dX = max_change(Xk, Xk1)
    dX = zeros(1,5);
    for i = 1:5, dX(i) = max(abs(Xk1{i}(:)-Xk{i}(:))); end
end

function init_e_values()
    sigma_fun = sigma_e_function();
    y0 = [r0;theta0;fai0;gamma0;polyfai0;V0];
    opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',de/2);
    [t,~] = ode45(@(e,y) dXde_ode(e,y,[0 0 0],sigma_fun), [e0,ef], y0, opts);
    idx = round(linspace(1,numel(t),N_points));
    e_values = t(idx); e_values = e_values(:);
    dei = zeros(N_points,1);  dei(2:end) = e_values(2:end)-e_values(1:end-1);
    dei_1 = zeros(N_points,1); dei_1(1:end-1) = e_values(2:end)-e_values(1:end-1);
end

function sigma_fun = sigma_e_function()
    rk=r0; thetak=theta0; faik=fai0; gammak=gamma0; polyfaik=polyfai0; Vk=V0; ek=e0; tk=0;
    dt = 1/sqrt(R0/g0); e_arr=[]; sig_arr=[];
    while ek <= ef
        sig = sigma_t(tk*sqrt(R0/g0)); cs=cosd(sig); ss=sind(sig);
        dx = dXdt(rk,thetak,faik,gammak,polyfaik,ek,cs,ss,0,0,0);
        rk=rk+dt*dx(1); thetak=thetak+dt*dx(2); faik=faik+dt*dx(3);
        gammak=gammak+dt*dx(4); polyfaik=polyfaik+dt*dx(5); Vk=Vk+dt*dx(6);
        ek=1/rk-Vk^2/2; tk=tk+dt; e_arr(end+1)=ek; sig_arr(end+1)=sig; %#ok<AGROW>
    end
    sigma_fun = @(e) interp1(e_arr, sig_arr, e, 'linear','extrap');
end

function s = sigma_t(t)
    if t<160, s=0; elseif t<600, s=-65; elseif t<1000, s=-50; else, s=45; end
end

function [rk,thetak,faik,gammak,polyfaik] = init_trajectories(sigma_fun)
    rk=zeros(N_points,Qcur); thetak=zeros(N_points,Qcur); faik=zeros(N_points,Qcur);
    gammak=zeros(N_points,Qcur); polyfaik=zeros(N_points,Qcur);
    opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',de);
    for k = 1:Qcur
        pert = Scur(k,:);
        y0 = [X0full(k,1);X0full(k,2);X0full(k,3);X0full(k,4);X0full(k,5);V0];
        [~,Y] = ode45(@(e,y) dXde_ode(e,y,pert,sigma_fun), e_values, y0, opts);
        rk(:,k)=Y(:,1); thetak(:,k)=Y(:,2); faik(:,k)=Y(:,3);
        gammak(:,k)=Y(:,4); polyfaik(:,k)=Y(:,5);
    end
end

function ydot = dXde_ode(e, y, pert, sigma_fun)
    sig = sigma_fun(e);
    ydot = dXde(y(1),y(2),y(3),y(4),y(5),y(6), cosd(sig), sind(sig), pert(1),pert(2),pert(3));
end

function ydot = dXde(r,theta,fai,gamma,polyfai,Vk,cs,ss,Rhodelta,CLdelta,CDdelta)
    e = 1/r - Vk^2/2;
    ydot = dXdt(r,theta,fai,gamma,polyfai,e,cs,ss,Rhodelta,CLdelta,CDdelta)/dedt(e,r,Rhodelta,CDdelta);
end

function de_val = dedt(e,r,Rhodelta,CDdelta)
    de_val = D_(e,r,CDdelta,Rhodelta)*V_(e,r);
end

function dx = dXdt(r,theta,fai,gamma,polyfai,e,cs,ss,Rhodelta,CLdelta,CDdelta)
    Vk=V_(e,r); Dk=D_(e,r,CDdelta,Rhodelta); Lk=L_(e,r,CLdelta,Rhodelta);
    cg=cos(gamma); sg=sin(gamma); cf=cos(fai); sf=sin(fai);
    cpf=cos(polyfai); spf=sin(polyfai); tf=tan(fai); tg=tan(gamma);
    dr=Vk*sg; dtheta=Vk*cg*spf/(r*cf); dfai=Vk*cg*cpf/r;
    dV=-Dk-sg/r^2+omega^2*r*cf*(sg*cf-cg*sf*cpf);
    dgamma=1/Vk*(Lk*cs+(Vk^2-1/r)*cg/r+2*omega*Vk*cf*spf+omega^2*r*cf*(cg*cf+sg*sf*cpf));
    dpolyfai=1/Vk*(Lk*ss/cg+Vk^2/r*cg*spf*tf+2*omega*Vk*(sf-tg*cpf*cf)+omega^2*r/cg*spf*sf*cf);
    dx=[dr;dtheta;dfai;dgamma;dpolyfai;dV];
end

function pc = path_constrain_init(Seta_cur, Qe)
    pc = zeros(N_points, Qe);
    for q = 1:Qe
        Rhodelta=Seta_cur(q,1); CLdelta=Seta_cur(q,2); CDdelta=Seta_cur(q,3);
        lQ = root_vec(@(r) Q_value_(e_values,r,Rhodelta), Qmax);
        lq = root_vec(@(r) q_value_(e_values,r,Rhodelta), qmax);
        ln = root_vec(@(r) n_value_(e_values,r,CLdelta,CDdelta,Rhodelta), gmax);
        pc(:,q) = max([lQ,lq,ln],[],2);
    end
end

function l = root_vec(fun, target)
    l = ones(N_points,1);
    f1 = fun(ones(N_points,1)); active = f1 > target;
    if ~any(active), return; end
    r_hi = min(r0, 1./e_values.*(1-1e-9));
    lo = ones(N_points,1); hi = r_hi;
    for iter = 1:80
        mid = 0.5*(lo+hi); fmid = fun(mid);
        sel1 = active & (fmid > target); sel2 = active & (fmid <= target);
        lo(sel1)=mid(sel1); hi(sel2)=mid(sel2);
    end
    l(active) = hi(active);
end

function [ok, r_n, t_n, f_n, g_n, p_n, u1n, u2n] = scp_iteration(rk, thetak, faik, gammak, polyfaik, pc, u1p, u2p)
    r = sdpvar(N_points,Qcur,'full'); theta = sdpvar(N_points,Qcur,'full');
    fai = sdpvar(N_points,Qcur,'full'); gamma = sdpvar(N_points,Qcur,'full');
    polyfai = sdpvar(N_points,Qcur,'full');
    cu1 = sdpvar(N_points,1,'full'); cu2 = sdpvar(N_points,1,'full');
    relax = sdpvar(2,1,'full');       % 终端区间均值到目标的松弛
    t_disp = sdpvar(2,1,'full');      % 终端区间方差上界（theta, phi）
    X = {r,theta,fai,gamma,polyfai}; u = {cu1,cu2};

    A = Ak(rk,thetak,faik,gammak,polyfaik);
    B = Bk(rk,thetak,faik,gammak,polyfaik);
    b = bk(rk,thetak,faik,gammak,polyfaik);

    Constraints = [];
    for i = 1:5
        Constraints = [Constraints, X{i}(1,:) == X0full(:,i).']; %#ok<AGROW>
    end
    Constraints = [Constraints, X{1}(end,:)*w_full == rf];
    % ---- 终端：标称均值到目标（软）+ 区间方差上界（AMM_4 Eq.30）----
    % 标称均值：在标称 eta=0 处对 xi 求期望（单点软约束，数值更稳）
    mu_theta = X{2}(end,:) * w_perf;    % 标量：标称 eta 处对 xi 的期望
    mu_fai   = X{3}(end,:) * w_perf;    % 标量
    Constraints = [Constraints, abs(mu_theta - thetaf) <= relax(1)];
    Constraints = [Constraints, abs(mu_fai   - faif)   <= relax(2)];
    Constraints = [Constraints, X{4}(end,:)*w_full >= gammaf(1), X{4}(end,:)*w_full <= gammaf(2)];
    Constraints = [Constraints, X{5}(end,:)*w_full >= polyfaif(1), X{5}(end,:)*w_full <= polyfaif(2)];
    % 区间方差上界：t_disp >= Var_xi[终端](eta_q) for all q（epigraph，凸）
    if Qcur > 1
        th_f = X{2}(end,:).';  ph_f = X{3}(end,:).';
        for q = 1:Qeta
            iq = (q-1)*Qxi + (1:Qxi);
            Constraints = [Constraints, t_disp(1) >= th_f(iq)'*M_xi*th_f(iq)]; %#ok<AGROW>
            Constraints = [Constraints, t_disp(2) >= ph_f(iq)'*M_xi*ph_f(iq)]; %#ok<AGROW>
        end
    end

    % 路径约束：对 xi 求期望（A_mu），对 eta 取最坏情况（每个 eta_q 都约束）
    r_mu = r * A_mu.';               % N x Q_eta
    rk_var = zeros(N_points, Qeta);
    for q = 1:Qeta
        idx = (q-1)*Qxi + (1:Qxi);
        rkq = rk(:,idx);
        rk_var(:,q) = sum(rkq .* (rkq*M_xi), 2);   % = diag(rkq*M_xi*rkq')
    end
    rk_var = max(rk_var, 0);
    Constraints = [Constraints, r_mu >= pc + kC*sqrt(rk_var)];

    for i = 1:5
        tmpi = zeros(N_points,Qcur); tmpi_1 = zeros(N_points,Qcur);
        for j = 1:5
            tmpi   = tmpi   + ((i==j) - A{i,j}.*dei/2).*X{j};
            tmpi_1 = tmpi_1 + ((i==j) + A{i,j}.*dei_1/2).*X{j};
        end
        for j = 1:2
            tmpi   = tmpi   - (B{i,j}.*dei/2).*(u{j}*ones(1,Qcur));
            tmpi_1 = tmpi_1 + (B{i,j}.*dei_1/2).*(u{j}*ones(1,Qcur));
        end
        tmpi   = tmpi   - b{i}.*dei/2;
        tmpi_1 = tmpi_1 + b{i}.*dei_1/2;
        Constraints = [Constraints, tmpi(2:end,:) == tmpi_1(1:end-1,:)]; %#ok<AGROW>
    end

    Constraints = [Constraints, cu1 >= cos(vmax), cu1 <= cos(vmin), cu1.^2 + cu2.^2 <= 1];
    Constraints = [Constraints, r >= 1, relax >= 0];
    Xk = {rk,thetak,faik,gammak,polyfaik};
    for i = 1:5
        Constraints = [Constraints, X{i}-Xk{i} <= delta(i)*tr_scale, X{i}-Xk{i} >= -delta(i)*tr_scale]; %#ok<AGROW>
    end

    Dk = D_(e_values,rk,Scur(:,3).',Scur(:,1).');
    Vk = V_(e_values,rk);
    Drk = Dr_(e_values,rk,Scur(:,3).',Scur(:,1).');
    objk = 1./(Dk.*Vk);  objk_r = -1./((Dk.*Vk).^2).*(Drk.*Vk);

    if Qcur > 1
        disp_obj = kJ * sum(t_disp);      % 终端区间方差上界（theta+phi）之和
    else
        disp_obj = 0;
    end
    % 性能项在标称 eta 处对 xi 求期望（w_perf），鲁棒项为区间方差上界 disp_obj
    obj = 0.5*trapz_e(polyfai*w_perf, e_values) ...
        + trapz_e((objk + objk_r.*(r - rk))*w_perf, e_values) ...
        + sum(relax) + disp_obj;

    options = sdpsettings('solver', solver_used, 'verbose',0,'debug',0, 'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 45);
    sol = optimize(Constraints, obj, options);
    if sol.problem ~= 0 && strcmpi(solver_used, 'mosek')
        % MOSEK 许可证失效/初始化失败时，回退 MATLAB 内置 coneprog
        sol = optimize(Constraints, obj, sdpsettings('solver','coneprog','verbose',0,'debug',0));
    end
    ok = (sol.problem == 0);
    if ~ok, warning('reentry_hybrid_mixed_pc:scp','SCO 子问题失败 problem=%d', sol.problem); end
    r_n=value(r); t_n=value(theta); f_n=value(fai); g_n=value(gamma); p_n=value(polyfai);
    u1n=value(cu1); u2n=value(cu2);
end

function val = trapz_e(f, e)
    val = 0.5*sum((f(1:end-1,:)+f(2:end,:)).*(e(2:end)-e(1:end-1)));
end

function A = Ak(r,theta,fai,gamma,polyfai)
    Drk = Dr_(e_values,r,Scur(:,3).',Scur(:,1).'); Dk = D_(e_values,r,Scur(:,3).',Scur(:,1).');
    sg=sin(gamma); cg=cos(gamma); sf=sin(fai); cf=cos(fai); tf=tan(fai);
    spf=sin(polyfai); cpf=cos(polyfai); Vk=V_(e_values,r);
    a11=-Drk.*sg./(Dk.^2); a14=cg./Dk;
    a21=-(Drk./(r.*Dk.^2)+1./(r.^2.*Dk)).*cg.*spf./cf;
    a23=cg.*spf.*sf./(r.*Dk.*cf.^2); a24=-sg.*spf./(r.*Dk.*cf); a25=cg.*cpf./(r.*Dk.*cf);
    a31=-(Drk./(r.*Dk.^2)+1./(r.^2.*Dk)).*cg.*cpf; a34=-sg.*cpf./(r.*Dk); a35=-cg.*spf./(r.*Dk);
    a41=Drk.*cg./(Dk.^2.*r).*(1./(r.*Vk.^2)-1)+cg./(Dk.*r.^2).*(2./(r.*Vk.^2)-1);
    a44=sg./(Dk.*r).*(1./(r.*Vk.^2)-1);
    a51=-cg.*spf.*tf.*(Drk./(r.*Dk.^2)+1./(r.^2.*Dk)); a53=cg.*spf./(r.*Dk.*cf.^2);
    a54=-sg.*spf.*tf./(r.*Dk); a55=cg.*cpf.*tf./(r.*Dk);
    z=zeros(size(a11));
    A={a11,z,z,a14,z; a21,z,a23,a24,a25; a31,z,z,a34,a35; a41,z,z,a44,z; a51,z,a53,a54,a55};
end

function B = Bk(r,theta,fai,gamma,polyfai)
    Dk=D_(e_values,r,Scur(:,3).',Scur(:,1).'); Lk=L_(e_values,r,Scur(:,2).',Scur(:,1).');
    cg=cos(gamma); Vk=V_(e_values,r);
    b41=Lk./Dk./(Vk.^2); b52=Lk./(Vk.^2.*cg.*Dk); z=zeros(size(b41));
    B={z,z; z,z; z,z; b41,z; z,b52};
end

function f0 = f0k(r,theta,fai,gamma,polyfai)
    Dk=D_(e_values,r,Scur(:,3).',Scur(:,1).');
    sg=sin(gamma); cg=cos(gamma); sf=sin(fai); cf=cos(fai); tf=tan(fai);
    spf=sin(polyfai); cpf=cos(polyfai); Vk=V_(e_values,r);
    f0={1./Dk.*sg; cg.*spf./(r.*Dk.*cf); cg.*cpf./(r.*Dk); (1./Dk-1./(r.*Vk.^2.*Dk)).*cg./r; cg.*spf.*tf./(r.*Dk)};
end

function fw = f_omega_k(r,theta,fai,gamma,polyfai)
    Dk=D_(e_values,r,Scur(:,3).',Scur(:,1).');
    sg=sin(gamma); cg=cos(gamma); sf=sin(fai); cf=cos(fai); tg=tan(gamma);
    spf=sin(polyfai); cpf=cos(polyfai); Vk=V_(e_values,r);
    b4=(2*omega*Vk.*cf.*spf+omega^2*r.*cf.*(cg.*cf+sg.*cpf.*sf))./(Vk.^2.*Dk);
    b5=(-2*omega*Vk.*(tg.*cpf.*cf-sf)+omega^2*r.*spf.*sf.*cf./cg)./(Vk.^2.*Dk);
    z=zeros(size(b4)); fw={z;z;z;b4;b5};
end

function b = bk(r,theta,fai,gamma,polyfai)
    f0=f0k(r,theta,fai,gamma,polyfai); fw=f_omega_k(r,theta,fai,gamma,polyfai); A=Ak(r,theta,fai,gamma,polyfai);
    x={r,theta,fai,gamma,polyfai}; b=cell(5,1);
    for i=1:5
        b{i}=f0{i}; tmp=zeros(size(f0{i}));
        for j=1:5, tmp=tmp+A{i,j}.*x{j}; end
        b{i}=b{i}-tmp+fw{i};
    end
end

function plot_all_trajectories(ROB_X, ROB_u)
    r=ROB_X{1}; theta=ROB_X{2}; fai=ROB_X{3}; gamma=ROB_X{4}; psi=ROB_X{5};
    Vv=V_(e_values,r);
    Dv=D_(e_values,r,Scur(:,3).',Scur(:,1).');
    t_e=cumtrapz(e_values,1./(Dv.*Vv),1);          % N x Q time (s)
    sig=atan2d(ROB_u{2},ROB_u{1});
    H=(r-1)*R0/1000;                              % altitude (km)
    Vkm=Vv*sqrt(g0*R0)/1000;                      % velocity (km/s)

    names_alt={'H (km)','\theta (deg)','\phi (deg)','\gamma (deg)','\psi (deg)','V (km/s)'};
    dat={H, theta*180/pi, fai*180/pi, gamma*180/pi, psi*180/pi, Vkm};

    figure('Name','Collocation states vs energy','Color','w');
    for k=1:6
        subplot(2,3,k); plot(e_values,dat{k},'LineWidth',0.3); xlabel('Energy e (non-dim)'); ylabel(names_alt{k}); grid on;
    end
    sgtitle('Collocation-point states versus energy');
    figure('Name','Control','Color','w');
    plot(e_values,sig,'LineWidth',1.5); xlabel('Energy e (non-dim)'); ylabel('Bank angle (deg)'); title('Bank-angle control (shared by all collocation points)'); grid on; drawnow;

    Qv=Q_value_(e_values,r,Scur(:,1).');
    qv=q_value_(e_values,r,Scur(:,1).');
    nv=n_value_(e_values,r,Scur(:,2).',Scur(:,3).',Scur(:,1).');
    figure('Name','Collocation path constraints','Color','w');
    subplot(3,1,1); hq=plot(e_values,Qv/1000,'LineWidth',0.3); hold on; hqm=yline(Qmax/1000,'k-.','LineWidth',1.5); ylabel('Heat flux\ $\dot{Q}_{\mathrm{max}}$\ (kW/m$^2$)','Interpreter','latex'); xlabel('Energy e (non-dim)'); legend([hq(1) hqm],{'trajectories','Q_{max}'},'Location','best'); grid on;
    subplot(3,1,2); hqp=plot(e_values,qv/1000,'LineWidth',0.3); hold on; hqpm=yline(qmax/1000,'k-.','LineWidth',1.5); ylabel('Dynamic pressure q (kPa)'); xlabel('Energy e (non-dim)'); legend([hqp(1) hqpm],{'trajectories','q_{max}'},'Location','best'); grid on;
    subplot(3,1,3); hn=plot(e_values,nv,'LineWidth',0.3); hold on; hnm=yline(gmax,'k-.','LineWidth',1.5); ylabel('Load factor n'); xlabel('Energy e (non-dim)'); legend([hn(1) hnm],{'trajectories','g_{max}'},'Location','best'); grid on;
    drawnow;
end

function [t_f,f_f,viol,traj] = mc_verify(u1c,u2c,N_mc)
    sigma_fun = @(e) atan2d(interp1(e_values,u2c,e,'linear','extrap'), interp1(e_values,u1c,e,'linear','extrap'));
    eta = 2*rand(N_mc,3)-1;   % eta：区间变量，在 [-1,1] 内均匀采样验证
    xi  = 2*rand(N_mc,2)-1;   % xi：均匀随机变量，在 [-1,1] 内均匀采样验证
    t_f=zeros(N_mc,1); f_f=zeros(N_mc,1); viol=zeros(N_mc,1);
    traj.r=zeros(N_points,N_mc); traj.theta=zeros(N_points,N_mc); traj.fai=zeros(N_points,N_mc);
    traj.gamma=zeros(N_points,N_mc); traj.psi=zeros(N_points,N_mc); traj.V=zeros(N_points,N_mc);
    traj.t=zeros(N_points,N_mc); traj.Q=zeros(N_points,N_mc); traj.q=zeros(N_points,N_mc); traj.n=zeros(N_points,N_mc);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',de);
    for k=1:N_mc
        pert = delta_eta .* eta(k,:);
        y0 = [r0; theta0+sigma_x0(1)*xi(k,1); fai0+sigma_x0(2)*xi(k,2); gamma0; polyfai0; V0];
        [~,Y] = ode45(@(e,y) dXde_ode(e,y,pert,sigma_fun), e_values, y0, opts);
        t_f(k)=Y(end,2); f_f(k)=Y(end,3); rk_mc=Y(:,1);
        traj.r(:,k)=Y(:,1); traj.theta(:,k)=Y(:,2); traj.fai(:,k)=Y(:,3);
        traj.gamma(:,k)=Y(:,4); traj.psi(:,k)=Y(:,5); traj.V(:,k)=Y(:,6);
        Vk=V_(e_values,rk_mc); Dk=D_(e_values,rk_mc,pert(3),pert(1));
        traj.t(:,k)=cumtrapz(e_values,1./(Dk.*Vk));
        traj.Q(:,k)=Q_value_(e_values,rk_mc,pert(1));
        traj.q(:,k)=q_value_(e_values,rk_mc,pert(1));
        traj.n(:,k)=n_value_(e_values,rk_mc,pert(2),pert(3),pert(1));
        viol(k)=any(traj.Q(:,k)>Qmax | traj.n(:,k)>gmax | traj.q(:,k)>qmax);
    end
end

function plot_mc_comparison(DO_X, DO_u, ROB_X, ROB_u, MC_DO, MC_ROB)
    % Nominal trajectory (ROB uses expectation over xi at eta=0)
    q_nom = find(all(abs(Seta)<1e-12,2),1);
    if isempty(q_nom), [~,q_nom]=min(sum(Seta.^2,2)); end
    idx_nom=(q_nom-1)*Qxi+(1:Qxi);
    r_DO=DO_X{1}(:,1); th_DO=DO_X{2}(:,1); fai_DO=DO_X{3}(:,1); gam_DO=DO_X{4}(:,1); psi_DO=DO_X{5}(:,1);
    r_ROB=ROB_X{1}(:,idx_nom)*w_xi; th_ROB=ROB_X{2}(:,idx_nom)*w_xi; fai_ROB=ROB_X{3}(:,idx_nom)*w_xi;
    gam_ROB=ROB_X{4}(:,idx_nom)*w_xi; psi_ROB=ROB_X{5}(:,idx_nom)*w_xi;
    H_DO=(r_DO-1)*R0/1000; H_ROB=(r_ROB-1)*R0/1000;
    V_DO=V_(e_values,r_DO)*sqrt(g0*R0)/1000; V_ROB=V_(e_values,r_ROB)*sqrt(g0*R0)/1000;

    % (1) Nominal optimal trajectories: 3D + altitude-velocity
    figure('Name','Nominal trajectory comparison','Color','w');
    subplot(1,2,1); plot3(th_DO*180/pi,fai_DO*180/pi,H_DO,'b-','LineWidth',2); hold on;
    plot3(th_ROB*180/pi,fai_ROB*180/pi,H_ROB,'r--','LineWidth',2);
    plot3(thetaf*180/pi,faif*180/pi,Hf/1000,'kp','MarkerSize',10,'MarkerFaceColor','k');
    xlabel('Longitude \theta (deg)'); ylabel('Latitude \phi (deg)'); zlabel('Altitude H (km)'); legend('DO','ROB','Target','Location','best'); grid on;
    subplot(1,2,2); plot(V_DO,H_DO,'b-','LineWidth',2); hold on; plot(V_ROB,H_ROB,'r--','LineWidth',2);
    set(gca,'XDir','reverse'); xlabel('Velocity (km/s)'); ylabel('Altitude (km)'); legend('DO','ROB','Location','best'); grid on;

    % (2) Monte Carlo 3D trajectories
    figure('Name','Monte Carlo 3D trajectories','Color','w');
    subplot(1,2,1); plot3(MC_DO.theta*180/pi, MC_DO.fai*180/pi, (MC_DO.r-1)*R0/1000,'LineWidth',0.3); grid on;
    xlabel('Longitude \theta (deg)'); ylabel('Latitude \phi (deg)'); zlabel('Altitude H (km)'); title('DO control');
    subplot(1,2,2); plot3(MC_ROB.theta*180/pi, MC_ROB.fai*180/pi, (MC_ROB.r-1)*R0/1000,'LineWidth',0.3); grid on;
    xlabel('Longitude \theta (deg)'); ylabel('Latitude \phi (deg)'); zlabel('Altitude H (km)'); title('ROB control');

    % (3) State curves: DO / ROB MC spread + nominal
    names={'Altitude H (km)','\theta (deg)','\phi (deg)','\gamma (deg)','\psi (deg)','V (km/s)'};
    dat_DO={(MC_DO.r-1)*R0/1000, MC_DO.theta*180/pi, MC_DO.fai*180/pi, MC_DO.gamma*180/pi, MC_DO.psi*180/pi, MC_DO.V*sqrt(g0*R0)/1000};
    dat_ROB={(MC_ROB.r-1)*R0/1000, MC_ROB.theta*180/pi, MC_ROB.fai*180/pi, MC_ROB.gamma*180/pi, MC_ROB.psi*180/pi, MC_ROB.V*sqrt(g0*R0)/1000};
    nom_DO={H_DO, th_DO*180/pi, fai_DO*180/pi, gam_DO*180/pi, psi_DO*180/pi, V_DO};
    nom_ROB={H_ROB, th_ROB*180/pi, fai_ROB*180/pi, gam_ROB*180/pi, psi_ROB*180/pi, V_ROB};
    figure('Name','State curves DO vs ROB','Color','w');
    for k=1:6
        subplot(2,3,k);
        h1=plot(e_values, dat_DO{k},'Color',[0.6 0.75 0.9],'LineWidth',0.4); hold on;
        h2=plot(e_values, dat_ROB{k},'Color',[0.95 0.7 0.5],'LineWidth',0.4);
        h3=plot(e_values, nom_DO{k},'b-','LineWidth',2);
        h4=plot(e_values, nom_ROB{k},'r--','LineWidth',2);
        xlabel('Energy e (non-dim)'); ylabel(names{k}); grid on;
        if k==1, legend([h1(1) h2(1) h3 h4],{'DO MC','ROB MC','DO nominal','ROB nominal'},'Location','best'); end
    end

    % (4) Path constraints: Q, q, n
    figure('Name','Path constraints DO vs ROB','Color','w');
    subplot(3,1,1);
    h1=plot(e_values, MC_DO.Q/1000,'Color',[0.6 0.75 0.9],'LineWidth',0.4); hold on;
    h2=plot(e_values, MC_ROB.Q/1000,'Color',[0.95 0.7 0.5],'LineWidth',0.4);
    h3=yline(Qmax/1000,'k-.','LineWidth',1.5); ylabel('Heat flux\ $\dot{Q}_{\mathrm{max}}$\ (kW/m$^2$)','Interpreter','latex'); xlabel('Energy e (non-dim)'); legend([h1(1) h2(1) h3],{'DO MC','RO MC','Q_{max}'},'Location','best'); grid on;
    subplot(3,1,2);
    h1=plot(e_values, MC_DO.q/1000,'Color',[0.6 0.75 0.9],'LineWidth',0.4); hold on;
    h2=plot(e_values, MC_ROB.q/1000,'Color',[0.95 0.7 0.5],'LineWidth',0.4);
    h3=yline(qmax/1000,'k-.','LineWidth',1.5); ylabel('Dynamic pressure q (kPa)'); xlabel('Energy e (non-dim)'); grid on;legend([h1(1) h2(1) h3],{'DO MC','RO MC','q_{max}'},'Location','best'); grid on;
    subplot(3,1,3);
    h1=plot(e_values, MC_DO.n,'Color',[0.6 0.75 0.9],'LineWidth',0.4); hold on;
    h2=plot(e_values, MC_ROB.n,'Color',[0.95 0.7 0.5],'LineWidth',0.4);
    h3=yline(gmax,'k-.','LineWidth',1.5); ylabel('Load factor n'); xlabel('Energy e (non-dim)'); grid on;legend([h1(1) h2(1) h3],{'DO MC','RO MC','g_{max}'},'Location','best'); grid on;

    % (5) Terminal dispersion
    figure('Name','Terminal dispersion','Color','w');
    scatter(MC_DO.theta(end,:)*180/pi, MC_DO.fai(end,:)*180/pi,18, 'b', 'filled', 'MarkerFaceAlpha',0.35); hold on;
    scatter(MC_ROB.theta(end,:)*180/pi, MC_ROB.fai(end,:)*180/pi,18, 'rs', 'MarkerFaceAlpha',0.35);
    plot(thetaf*180/pi, faif*180/pi,'kp','MarkerSize',18,'MarkerFaceColor','k');
    xlabel('Longitude \theta (deg)'); ylabel('Latitude \phi (deg)'); legend('DO','ROB','Target','Location','best'); grid on; drawnow;
end

%% ===================== 求积与气动/动力学模型 =====================

function [x,w] = gauss_legendre(n)
    beta=(1:n-1)./sqrt(4*(1:n-1).^2-1); J=diag(beta,1)+diag(beta,-1);
    [V,D]=eig(J); x=diag(D).'; w=2*V(1,:).^2; [x,o]=sort(x); w=w(o); w=w/sum(w);
end

function [x,w] = cheb_lobatto(n)
    k=0:n-1; x=cos(k*pi/(n-1)); w=ones(1,n); w(1)=0.5; w(end)=0.5; w=w/sum(w);
end

function [pts,w] = tensor_quad(nodes_cell, weights_cell)
    d=numel(nodes_cell); [G{1:d}]=ndgrid(nodes_cell{:}); [W{1:d}]=ndgrid(weights_cell{:});
    Q=prod(cellfun(@numel,nodes_cell)); pts=zeros(Q,d); w=ones(Q,1);
    for i=1:d, pts(:,i)=G{i}(:); w=w.*W{i}(:); end
end

function Vk = V_(e,r), Vk=sqrt(2*(1./r-e)); end
function rho = Rho_(r,Rhodelta), rho=1.225*exp(-(r-1)*R0*13.785e-5).*(1+Rhodelta); end
function rho_r = Rho_r_(r,Rhodelta), rho_r=-Rho_(r,Rhodelta)*R0*13.785e-5; end
function a = sound_speed_(r)
    T=calculate_temperature((r-1)*R0/1000); a=(331.6+0.6*(T-273.15))/sqrt(g0*R0);
end
function T = calculate_temperature(h)
    T=288.15*ones(size(h)); layers=[0,-6.5;11,0;20,1;32,2.8;47,0;51,-2.8;71,-2;85,NaN];
    for i=1:size(layers,1)-1
        m1=(h>=layers(i,1))&(h<layers(i+1,1)); m2=(h>=layers(i,1))&(h>=layers(i+1,1));
        T=T+m1.*(layers(i,2).*(h-layers(i,1)))+m2.*(layers(i,2).*(layers(i+1,1)-layers(i,1)));
    end
end
function alpha = alpha_profile(Mach), alpha=45-0.612*max(10-Mach,0).^2; end
function CL = CL_(alpha,mach,CLdelta), CL=Calcu_CLCD(alpha,mach,false).*(1+CLdelta); end
function CD = CD_(alpha,mach,CDdelta), CD=Calcu_CLCD(alpha,mach,true).*(1+CDdelta); end

function out = Calcu_CLCD(alpha,mach,isCD)
    mask=mach<=5;
    if isCD
        out=zeros(size(alpha));
        a=0.1908219490086461; b=0.7921786756723015; c=81.77767906793568; d=0.7608546899234436; e=-0.06715377897523668; f=4.809129538309562; g=-0.8220180925372075;
        z=a+b*exp(-0.5*(log(alpha/c)/d).^2)+e*exp(-0.5*((mach-f)/g).^2); out(mask)=z(mask);
        a=0.2893865934644147; b=0.7101395427838309; c=44.69254972733151; d=12.43577575552023; e=-0.243653871914911; f=5.428869369818743; g=0.2330424954375339; hh=1.098649201604347;
        z1=1./(1+exp(-(alpha-c)/d)); z2=1./(1+exp(-(mach-f)/g)); z=a+b*z1+e*z2+hh*z1.*z2; out(~mask)=z(~mask);
    else
        out=zeros(size(alpha));
        a=0.6736442631404929; b=-0.2396729544298644; c=0.02354634960855771; d=0.5137518146743224; e=35.7738334917245; f=0.5482342780270542;
        z1=exp(-0.5*(log(alpha/e)/f).^2); z=a+mach.*(b+c*mach)+d*z1; out(mask)=z(mask);
        a=-0.2361725195924329; b=1.037601261644054; c=50.11901193898274; d=26.17064737989058; e=0.8761093802115429; f=4.857160518263465; g=-0.485173805363233; hh=-1.164121063690802;
        z1=exp(-0.5*((alpha-c)/d).^2); z2=exp(-0.5*((mach-f)/g).^2); z=a+b*z1+e*z2+hh*z1.*z2; out(~mask)=z(~mask);
    end
end

function L = L_(e,r,CLdelta,Rhodelta)
    Mach=V_(e,r)./sound_speed_(r); L=0.5*R0*Rho_(r,Rhodelta).*V_(e,r).^2*Sref.*CL_(alpha_profile(Mach),Mach,CLdelta)/m;
end
function D = D_(e,r,CDdelta,Rhodelta)
    Mach=V_(e,r)./sound_speed_(r); D=0.5*R0*Rho_(r,Rhodelta).*V_(e,r).^2*Sref.*CD_(alpha_profile(Mach),Mach,CDdelta)/m;
end
function Dr = Dr_(e,r,CDdelta,Rhodelta)
    Mach=V_(e,r)./sound_speed_(r);
    dsdr=num_deriv(@(rr) sound_speed_(rr),r,1e-7); dMach=-V_(e,r)./(sound_speed_(r).^2).*dsdr;
    dCD=num_deriv(@(mm) CD_(alpha_profile(mm),mm,CDdelta),Mach,1e-7);
    Dr=0.5*R0*Rho_r_(r,Rhodelta).*V_(e,r).^2*Sref.*CD_(alpha_profile(Mach),Mach,CDdelta)/m ...
      +0.5*R0*Rho_(r,Rhodelta).*V_(e,r).^2*Sref.*dCD.*dMach/m;
end
function d = num_deriv(fun,x,dx), d=(fun(x+dx)-fun(x-dx))/(2*dx); end
function n = n_value_(e,r,CLdelta,CDdelta,Rhodelta), n=sqrt(L_(e,r,CLdelta,Rhodelta).^2+D_(e,r,CDdelta,Rhodelta).^2); end
function q = q_value_(e,r,Rhodelta), q=0.5*g0*R0*Rho_(r,Rhodelta).*V_(e,r).^2; end
function Q = Q_value_(e,r,Rhodelta), Q=kQ*sqrt(g0*R0)^3.15.*sqrt(Rho_(r,Rhodelta)).*V_(e,r).^3.15; end

function P = legpoly(k, x)
% legpoly  归一化 Legendre 多项式（U(-1,1)），满足 E[P_k P_l] = delta_kl。
    if k == 0, P = ones(size(x)); return; end
    if k == 1, P = sqrt(3)*x; return; end
    P0 = ones(size(x)); P1 = x;
    for i = 2:k
        P2 = ((2*i-1)*x.*P1 - (i-1)*P0)/i; P0 = P1; P1 = P2;
    end
    P = sqrt(2*k+1)*P1;
end


end
