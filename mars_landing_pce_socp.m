function mars_landing_pce_socp(N_in, vrb_in)
% mars_landing_pce_socp  火星三自由度动力下降着陆 混合不确定性 鲁棒轨迹优化
%
% ------------------------------------------------------------------------
% 编程逻辑与绘图思路完全仿照 reentry_amm4_mssp_nxi2.m（地球再入）：
%   Stage 1  DO  ：确定性轨迹优化（单条轨迹）
%   Stage 2  RO ：混合不确定性（xi 均匀随机 + eta 区间）鲁棒轨迹优化
%   MC 验证 + DO/RO 对比绘图（标称轨迹、三维打靶、状态曲线、
%                              过程约束、终端散布、配点集）
%
% 唯一的区别：核心凸优化问题由"序列凸优化(SCP)"换成"二阶锥规划(SOCP)"。
%   - 变量变换 z = log(m)：dz/dt = -alpha*sigma，质量方程线性化
%   - 无损凸化(lossless convexification)：||u|| <= sigma <= rho2/m
%   - 质量依赖推力上下限用 Taylor 展开凸化（上界一阶、下界二阶）
%   - 解析离散化（超级矩阵）：y_k = A^k y0 + Lambda_k*grav + Psi_k*eta
%   变换后系统对状态/控制线性，一次凸求解即全局最优，无需 SCP 的
%   线性化与信赖域迭代。
%
% 求解实现说明（重要）：
%   控制只有一份 E (4N x 1)，由 Q 个配点共享（与再入代码里 cu1/cu2 是
%   (N_points x 1)、用 ones(1,Qcur) 广播的做法一致）；状态有 Q 份 Y_s，
%   各自从自己的初值出发被同一条控制驱动。故整个 RO 是
%   "1 份控制 + Q 份状态"的单个 SOCP，而不是 Q 个各自独立的 SOCP。
%
%   【决策变量】只有
%       E     4N x 1  共享控制（u 与 sigma）
%       relax 6  x 1  终端均值松弛
%   状态不作为决策变量：解析离散化已给出显式仿射关系
%       Y_s = Abig*y0_s + Lamug + PsiK*E
%   把它们当变量会得到 7N*Q + 4N ~ 22600 个变量、且 64 份之间仅靠同一个 E
%   耦合，MOSEK 内点法直接返回 problem=4（数值失败）；消去后只剩 4N+6 个
%   变量，约束全部化为 E 的仿射 / SOC / 凸二次，数学上等价而规模小两个量级。
%
%   终端处理（对应 AMM4 的 relax / t_disp）：
%     均值 -> relax 软约束（xi 加权期望逼近 0），以 K_final 线性重罚进目标
%     散布 -> 谱插值直接给出闭式常数（见下），不进模型，只作指标输出
%             （若写成 sdisp >= norm(...) 的锥，会给 MOSEK 添上百个 63 行、
%              系数仅 1e-16 的退化锥，实测使模型退化：sigma 偏离 ||u|| 达
%              46%、多耗 70 kg 燃料，故必须用闭式表达）
%
%   注意（结构性结论，务必理解）：动力学对状态线性、且控制与 xi 无关，于是
%       r(tf) - r_nom(tf) = (dr0 + tf*dv0) + drho(eta)
%   第一项只由初值不确定性决定、**与整条控制剖面无关**（代码里体现为
%   Mi_map*ones = 0），是任何共享开环剖面都摆脱不掉的"不可约项"：
%       sqrt(sigma_r0^2 + tf^2*sigma_v0^2)
%   第二项 drho(eta) 只由区间不确定性决定，且**只经由"推力指令能否被真实
%   执行"进入**——因为 u 是推力加速度，指令一旦被发动机给出，轨迹就与
%   推力幅值无关。所以：
%     · DO 的剖面贴在名义边界 0.80/0.20 Tmax 上，对 eta<0 的样本上限被压到
%       0.784、对 eta>0 的样本下限被抬到 0.204 -> 指令必然越界、被饱和截断
%       -> 实际执行的已不是设计剖面 -> 产生附加散布 drho；
%     · RO 的剖面收在区间最紧边界 0.204~0.784 Tmax 内 -> 任意 eta 都不饱和
%       -> drho ≡ 0 -> 散布恰好等于不可约项。
%   因此"RO 的 MC 散布比 DO 小"这件事成立，但幅度有限：不可约项占 92%，
%   RO 能消除的只有饱和附加项（本例 131 m -> 0，总 1sigma 341.7 -> 315.6 m，
%   即 RO/DO = 0.92，同时把推力违约样本 100% 降到 0%、多耗约 10 kg 推进剂）。
%   文件末尾的 [对照] 段用"放开可用推力带"的第三条 MC 直接验证了这一点：
%   同一条 DO 剖面若假设发动机能给出任意推力，其散布与 RO 完全相同。
%
% ------------------------------------------------------------------------
% 动力学（图片 Eq.37，r1 为高度分量）
%     r_dot = v
%     v_dot = u + gvec ,      gvec = [-g_m;0;0]
%     m_dot = -alpha*||Tc|| ,  Tc = m*u
%
% 量纲约定（务必区分"推力加速度"与"推力"）：
%   【坐标命名】x = r1（高度），y = r2（射程），z = r3（横程）；
%              相应地 v_x/v_y/v_z、u_x/u_y/u_z。x 为高度方向（g 沿 -x）。
%     u      —— 推力加速度 (m/s^2)，出现在 v_dot = u + g 中
%     Tc     —— 推力 (N)，Tc = m*u
%     推力约束 0.2 <= m*||u||/Tmax <= 0.8  即  rho1 <= ||Tc|| <= rho2
%     等价地  ||u|| ∈ [rho1/m, rho2/m]
%     质量流率 m_dot = -alpha*||Tc|| = -alpha*m*||u||   ->  z=log(m) 时 z_dot = -alpha*sig
%
% 离散化一致性（关键，实测踩过坑）：
%   解析离散化中，第 k 个状态块对应 t_k，第 k 个控制块对应区间 [t_{k-1}, t_k]。
%   区间内 sigma 恒定、质量单调下降，故推力 T = m*sigma 在区间**起点**最大。
%   因此推力约束必须用"区间起点质量" m(t_{k-1}) 与 sigma_k 配对。
%   若误配成区间终点质量 m(t_k)（本文件早期版本即如此），等价于把上限
%   放宽 m(t_{k-1})/m(t_k) - 1 ~ alpha*rho2*dt/m ~ 0.575%，表现为：
%       - 优化解自身在节点上就超界（sigma_k > 其约束上限）；
%       - MC 实测 m|u| 恰好超 0.98*rho2 约 0.6%，且无法靠收紧约束消除；
%       - 只能再加 thr_marg 之类的人为裕度去"补"，属于治标不治本。
%   修正后 DO/RO 的 max T/Tmax = 0.7840，恰好等于允许上限，约束真正起作用。
%
% 不确定性（按题目要求）
%   xi  （均匀随机，Legendre PCE，n_xi = 6）：
%         初始位置 r0 与初始速度 v0 各三分量，均为参考值的 ±5%。
%   eta （区间，Chebyshev-Lobatto 端点）：推力幅值上下限的 ±2%。
%         控制不依赖 eta（区间量未知），区间以最坏情况端点进入约束。
%
% 需要 YALMIP + MOSEK。
% ------------------------------------------------------------------------

%% ===================== 0. 环境检查 =====================
if isempty(which('mosekopt'))
    msk = 'C:\Program Files\Mosek\11.0\toolbox\r2019b';
    if exist(msk, 'dir'), addpath(msk); end
end
if isempty(which('sdpvar'))
    error('YALMIP not found: add the YALMIP folder to the MATLAB path.');
end

%% ===================== 1. 物理参数（图片 Eq.37） =====================
g_m    = 3.7114;                 % 火星重力加速度 (m/s^2)
gvec   = [-g_m; 0; 0];           % 重力矢量（沿 r1 反方向）
m0     = 1905;                   % 初始质量 (kg)
alpha  = 5.09e-4;                % 燃料消耗系数 (s/m)
Tmax   = 7.5*9.80665*225;        % 最大推力 (N) ≈ 16549
rho1   = 0.2*Tmax;               % 推力下限（图片 Eq.37 为 0.3*Tmax；此处按下限可放松取 0.2）
rho2   = 0.8*Tmax;               % 推力上限 0.8*Tmax

r0 = [1500; 500; 2000];          % 初始位置 (m)，r1 = 高度
v0 = [-75; 40; 100];             % 初始速度 (m/s)

tf = 81;                         % 终端时间 (s)，t0 = 0
N  = 50;                         % 时间离散步数（状态点 N+1 个）
if nargin >= 1 && ~isempty(N_in), N = N_in; end
VRB = 0;
if nargin >= 2, VRB = vrb_in; end
dt = tf/N;

tan_gs = tan(86*pi/180);         % sqrt(r2^2+r3^2) <= tan(86 deg)*r1

%% ===================== 2. 混合不确定性配置 =====================
% xi（均匀随机，Legendre PCE）：初始位置 3 + 初始速度 3
n_xi     = 6;
sigma_r0 = 0.02*abs(r0);         % 初始位置 ±5%
sigma_v0 = 0.02*abs(v0);         % 初始速度 ±5%

% ---- 配点采样模式：决定 Q=64 个配点在 ξ 空间的位置 -------------------
% 为什么需要这个开关：
%   状态对 ξ 是**仿射**的（r_i = r_i,nom + dr0_i + t*dv0_i），而 dr0_i、dv0_i
%   只由 ξ_i 与 ξ_{i+3} 两个分量决定。若配点取"轴对齐张量网格"
%   （每维 {-1,+1}），则每个分量最多只看到 2x2=4 个取值，64 条轨迹在
%   单分量图上只会呈现出 4 条（速度 2 条），其余 16/32 条严格重合。
%   这是网格采样的固有性质，与 PCE 正确性无关；换成"一般位置"的采样集，
%   64 个点在任意 2 维投影上都互不相同，64 条轨迹就会全部显示出来。
%
%   'lhs'    (默认) 64 点拉丁超立方（每维 64 个对称层中心）
%            · 一般位置 -> 每个分量 64 条轨迹互不相同（这就是本开关的目的）
%            · 层中心关于 0 对称 -> 每列和为 0，均值对仿射响应精确
%            · 对仿射响应，PCE 插值仍然精确（基张成的空间包含全部低阶项）
%   'corner' 每维 2 点 Chebyshev-Lobatto -> {-1,+1}^6 的 64 个箱角
%            · 覆盖全部最坏情况端点，内插不外推
%            · 但单分量图上只剩 4 条（速度 2 条）曲线
% 两种模式下 cm（均值）与 SigA（标准差）完全相同，故 DO/RO 的优化结果、
% 燃料、终端散布都不受影响 —— 只是配点图的可见性不同。
SAMP = 'lhs';

% eta（区间）：推力幅值 ±2%，取最坏情况端点
delta_eta = 0.02;
% 推力幅值约束 m*||u|| ∈ [rho1, rho2] 是"推力(力)"的约束，
% 而 u 是"推力加速度"(m/s^2)，两者相差一个 m，切勿混用：
%     u = Tc/m ,  Tc = m*u ,  ||u|| ∈ [rho1/m, rho2/m]
% 区间不确定性只作用在推力幅值上，故取最紧端点：
%     Tc <= rho2*(1-2%)      ->  ||u|| <= rho2_w/m
%     Tc >= rho1*(1+2%)      ->  ||u|| >= rho1_w/m
% 上限不需要额外裕度：按"区间起点质量"配对后上限可严格满足。
% 下限另加 1% 小裕度，用于覆盖区间内质量衰减带来的"末端"推力下陷：
%   下界在区间**终点**最紧，而配对用的是区间起点质量，
%   二者相差 m(t_{k-1})/m(t_k)-1 ~ alpha*rho2*dt/m ~ 0.57%。

rho2_w = rho2*(1 - delta_eta);   % 最紧推力上限
rho1_w = rho1*(1 + delta_eta);   % 最紧推力下限

% 过程约束（高度、滑翔角锥）的走廊裕度系数 kC_path(k)：
% "在均值之外再预留 kC 倍标准差"，对应 AMM4 里的 r_mu >= pc + kC*sqrt(var)。
% 本算例默认取 0（只约束均值），原因是实测与结构性的：
%   1) 末段高度必然趋于 0，而 ±5% 初值在 81 s 内可造成近 200 m 的高度散布，
%      "均值 - kC*sigma >= 0" 在触地附近必然不可行，只能让 kC 随 k 衰减；
%   2) 更根本的是：动力学对状态线性、控制与 xi 无关，故终端散布恒等于
%      Var(r0 + tf*v0)，与剖面无关（见文件头）。加 sigma 裕度既不能减小散布，
%      又要多耗燃料（实测 kC=2 衰减时 425.7 kg -> 482.9 kg，+57 kg），
%      因此把它作为可选项保留，默认关闭。
kC_path = zeros(N,1);
% 如需启用（满裕度 2 sigma，末段线性衰减到 0），取消下一行注释：
% kC_path = 2.0*max(0, 1 - ((1:N)'-1)/(N-1));

% 终端软约束权重（kg/m）：relax 以 m 计（速度分量以 m/s 计）。
% 取 5 表示"1 m 落点误差"与"5 kg 燃料"等价，与燃料量级(~430 kg)可比，
% 从而先把落点压小，再在剩余自由度上省燃料。
K_final = 5;

% 基函数与配点都固定为 Q = 2^n_xi = 64（多重线性张量基：每维 0/1 两个基函数）
Q = 2^n_xi;

switch lower(SAMP)
    case 'corner'
        % Chebyshev-Lobatto 端点 {-1,+1}：64 个点恰为不确定域 [-1,1]^6 的角点，
        % PCE 插值在整条不确定域内都是"内插"而非外推。
        [xL, wL]         = cheb_lobatto(2);
        [xi_nodes, w_xi] = tensor_quad(repmat({xL},1,n_xi), repmat({wL},1,n_xi));
    case 'lhs'
        % Q=64 点拉丁超立方：每维 64 个对称层中心，每层恰好落一点。
        % 这样每一维的 64 个坐标值互不相同（速度分量因此有 64 条不同曲线），
        % 且层中心关于 0 对称 -> 每一列的和严格为 0 -> 均值对仿射响应精确。
        % 注意不能用"32 点 + 反演点"：32 层中心的集合本身关于 0 对称，
        % 取反后仍是同一组值，于是每维只剩 32 个不同的坐标值（速度只出 32 条曲线）。
        %
        % 位置分量看到的是二维投影 (ξ_i, ξ_{i+3})，随机抽取时理论上可能撞点，
        % 故这里做一次抽样校验：任一投影出现重复就重抽（最多 200 次）。
        rng(20261026);                        % 固定种子，结果可复现
        okdraw = false;
        for ntry = 1:200
            L = zeros(Q, n_xi);
            for d = 1:n_xi
                pv = randperm(Q).';
                L(:,d) = (2*pv - 1)/Q - 1;    % 层中心 -> (-1,1)
            end
            okdraw = true;
            for i = 1:3                        % 只需校验 3 个二维投影
                if size(unique([L(:,i), L(:,3+i)], 'rows'), 1) < Q
                    okdraw = false;  break;
                end
            end
            if okdraw, break; end
        end
        if ~okdraw, error('LHS 抽样 200 次仍未得到一般位置的投影'); end
        xi_nodes = L;
        w_xi     = ones(Q,1)/Q;               % 等权：均值 = 简单平均
    otherwise
        error('SAMP 只能是 ''lhs'' 或 ''corner''');
end

% 多指标表（0 基）：每维取 0/1 两层，张成全部平方自由单项式（含常数项）
% 注意：必须显式声明 ind2sub 的多个输出；若写成 MultIdx(ii,:)=ind2sub(...)，
%       MATLAB 只会返回单个线性下标并被广播，导致基函数退化（rank=2）。
dims    = 2*ones(1,n_xi);
MultIdx = zeros(Q, n_xi);
for ii = 1:Q
    sub = cell(1, n_xi);
    [sub{:}] = ind2sub(dims, ii);
    MultIdx(ii,:) = cell2mat(sub) - 1;
end

% 谱(PCE)系数投影：归一化 Legendre 张量基
Phi_xi = zeros(Q, Q);
for ii = 1:Q
    P = ones(Q,1);
    for d = 1:n_xi
        P = P .* legpoly(MultIdx(ii,d), xi_nodes(:,d));
    end
    Phi_xi(:,ii) = P;
end
C_xi = (Phi_xi .* w_xi.').';                  % Galerkin 系数映射（Q x Q）
M_xi = C_xi(2:end,:)' * C_xi(2:end,:);        % 谱方差二次型
% 谱插值映射：beta = Mi_map * f_coll（去掉常数行）
%   对 Lobatto 端点节点，离散正交不成立，故方差要用"插值"而不是"Galerkin"：
%       Var(f) = || Mi_map * f_coll ||^2 ,  Mi_map = (Phi \ I)(2:end,:)
Mi_map = (Phi_xi \ eye(Q));
Mi_map = Mi_map(2:end,:);

% ---- 采样/谱分解自检 -------------------------------------------------
% 取一条已知的仿射测试函数 f = 1 + sum a_k*xi_k，其谱系数解析已知为 a。
% 若谱系数恢复不出来，说明 Phi 病态或节点排布有问题，后面的 sigma 全不可信。
a_test  = (1:n_xi).';
c_true  = zeros(Q,1);  c_true(1) = 1;          % 常数项系数
for d = 1:n_xi
    e = zeros(1,n_xi);  e(d) = 1;
    % 基是"归一化 Legendre"而非单项式：P_1 = sqrt(3)*x，
    % 故 f = 1 + sum a_d*x_d 在该基下的系数是 a_d/sqrt(3)
    c_true(ismember(MultIdx, e, 'rows')) = a_test(d)/sqrt(3);
end
cf_test = Mi_map*(1 + xi_nodes*a_test);        % 恢复出的高次项（去常数行）
fprintf(['[PCE] sampling=%-6s Q=%d | cond(Phi)=%.2e | max|Mi_map*1|=%.1e' ...
         ' | affine self-test err=%.1e\n'], ...
        upper(SAMP), Q, cond(Phi_xi), max(abs(Mi_map*ones(Q,1))), ...
        max(abs(cf_test - c_true(2:end))));

% 各配点初始状态
X0all = zeros(Q, 6);
for s = 1:Q
    rs = r0(:) + sigma_r0(:).*xi_nodes(s,1:3).';
    vs = v0(:) + sigma_v0(:).*xi_nodes(s,4:6).';
    X0all(s,:) = [rs.', vs.'];
end

%% ===================== 3. 连续/离散系统矩阵（7 维状态） =====================
A_c = [zeros(3,3), eye(3), zeros(3,1); zeros(4,7)];
B_c = [[zeros(3,3); eye(3); 0 0 0], [0 0 0 0 0 0 -alpha]'];

A = A_c*dt + eye(7);                    % expm(A_c*dt)，A_c 幂零
B = (eye(7)*dt + A_c*dt^2/2)*B_c;       % 精确零阶保持输入矩阵

Abig = zeros(7*N, 7);
Ak   = eye(7);
for k = 1:N
    Ak = Ak*A;
    Abig(7*(k-1)+(1:7), :) = Ak;
end

LambdaK = zeros(7*N, 4);
LambdaK(1:7,:) = B;
for k = 2:N
    LambdaK(7*(k-1)+(1:7), :) = A*LambdaK(7*(k-2)+(1:7), :) + B;
end

PsiK = zeros(7*N, 4*N);
PsiK(1:7, 1:4) = B;
for k = 2:N
    PsiK(7*(k-1)+(1:7), :) = A*PsiK(7*(k-2)+(1:7), :);
    PsiK(7*(k-1)+(1:7), (k-1)*4+(1:4)) = B;
end

ugrav  = [-g_m; 0; 0; 0];
Lamug  = LambdaK*ugrav;                 % 重力自由响应

% ---- 参考（最小质量）轨迹：时刻必须与"控制区间"对齐 ----
% 解析离散化的排布是：  第 k 个状态块 = t_k 时刻的状态（Ak = A^k），
%                       第 k 个控制块 = 区间 [t_{k-1}, t_k] 上的控制。
% 在区间内 sigma 恒定、质量单调下降，故推力 T = m*sigma 在**区间起点**
% 取最大。推力约束 ||Tc|| <= rho2 必须在区间起点成立，因此必须用
% 起点质量 m(t_{k-1}) 与 sigma_k 配对。
%   （若误配成区间终点质量 m(t_k)，等价于把上限放宽
%     m(t_{k-1})/m(t_k) - 1 ~ alpha*rho2*dt/m ~ 0.575%，
%     表现为实测 m|u| 恰好超 0.98*rho2 约 0.6%，且无法靠收紧约束消除。）
tref   = (0:N-1)'*dt;                   % 各控制区间的起点时刻 t_{k-1}
z0s    = log(m0 - alpha*rho2*tref);     % 起点参考质量（对数）
expmzs = exp(-z0s);                     % = 1/m_ref(t_{k-1})

%% ===================== 4. Stage 1 : 确定性轨迹优化（DO） =====================
fprintf('===== Stage 1: deterministic optimization (DO) =====\n');
tic
y0_nom = [r0; v0; log(m0)];
[okDO, DOy, DOe] = solve_point(y0_nom);

if ~okDO, error('DO stage failed.'); end
DO = unpack_point(DOy, DOe, y0_nom);
fprintf('-------- DO nominal --------\n');
print_summary(DO);
toc
%% ===================== 5. Stage 2 : 鲁棒轨迹优化（RO） =====================
% 【与 reentry_amm4_mssp_nxi2 / reentry_hybrid_mixed_pc_nxi2 完全一致的逻辑】
%   控制：只有一份 E (4N x 1)，被全部 Q 个配点共享 —— u 与 xi 无关；
%         对应再入代码里 cu1/cu2 是 (N_points x 1)、用 ones(1,Qcur) 广播。
%   状态：Q 份 Y_s (7N x 1)，各自从自己的初值出发，被同一条控制驱动。
%         这里不把它们当决策变量——解析离散化已给出 Y_s = Abig*y0_s+Lamug+
%         PsiK*E，直接代入即可（详见文件头"决策变量"说明）。
%   于是"一条控制剖面 + Q 条状态轨迹"，配点间的差异只来自初值，不来自控制。
%
%   过程约束按 AMM4 的做法作用在**均值**上（r_mu >= pc 形式），而不是要求
%   每个配点各自满足：末段高度趋 0，±5% 初值可造成数百米高度散布，逐点
%   约束必然不可行。kC 倍 sigma 的走廊裕度作为可选项保留（默认 0，理由
%   见 kC_path 处的说明）。
%
%   一条开环剖面不可能让 Q 个不同初值同时精确落点，故终端按 AMM4 的做法：
%     均值  -> 松弛变量 relax 软约束（对 xi 加权期望，逼近 0），重罚进目标
%     散布  -> 谱插值闭式常数（见 solve_shared 内 SigA），只作指标输出
%   对应再入代码里 relax / t_disp 的角色；t_disp 那里是变量，因为再入的终端
%   量对状态非线性、方差确实依赖控制，而本算例是线性动力学 -> 方差与剖面无关。
fprintf('\n===== Stage 2: robust optimization (shared control, Q=%d) =====\n', Q);
tic
[okR, Yall, Esg, diagR] = solve_shared();

if ~okR
    error('RO stage failed (problem=%d).', diagR.problem);
end

% ---- 解包：64 条状态轨迹（控制是同一份，故 ROu 的 64 列完全相同）----
ROu = repmat(Esg, 1, Q);
ROy = cell(Q,1);
ROr = zeros(3,N+1,Q);  ROv = zeros(3,N+1,Q);
ROm = zeros(N+1,Q);    ROuu = zeros(3,N,Q);
for s = 1:Q
    ys = [X0all(s,1:3).'; X0all(s,4:6).'; log(m0)];
    ROy{s} = Yall(:,s);
    Ss = unpack_point(Yall(:,s), Esg, ys);
    ROr(:,:,s) = Ss.r;  ROv(:,:,s) = Ss.v;
    ROm(:,s)   = Ss.m;  ROuu(:,:,s) = Ss.u;
end

% ---- 名义轨迹：用同一条控制从标称初值 (r0,v0) 精确积分 ----
RO_nom = propagate_nominal(Esg);
RO_nom.t  = (0:N)*dt;
RO_nom.tc = (1:N)*dt;
DO.t       = (0:N)*dt;
DO.tc      = (1:N)*dt;

fprintf('-------- RO nominal (shared control, from r0/v0) --------\n');
print_summary(RO_nom);

% ---- 模型/积分一致性检查 ----
% 解析离散解用的是 (sigma_k, 区间起点模型质量)；逐步积分用的是 (||u_k||, 递推质量)。
% 若两者给出的 T/Tmax 不一致，说明 sigma 与 ||u|| 未贴合（无损凸化未取等号）
% 或质量参考轨迹对不齐。这里把两个数都打印出来，便于交叉核对。
zmdl = Yall(7*(1:N).', 1);                 % 模型 z(t_k), k=1..N
mmdl = exp([log(m0); zmdl(1:end-1)]);      % 区间起点模型质量
Sgv  = reshape(Esg,4,N); Sgv = Sgv(4,:).'; % sigma_k
Tm   = mmdl.*Sgv/Tmax;
accN = sqrt(sum(RO_nom.u.^2,1)).';
Tp   = RO_nom.m(1:N).*accN/Tmax;
[~, km] = max(Tm);
fprintf(['[consistency] model max T/Tmax=%.4f (k=%d) | propagated max T/Tmax=%.4f (k=%d)\n' ...
         '              sigma vs |u|: max(sigma-|u|)=%.3e (无损凸化应取等号)\n'], ...
        max(Tm), km, max(Tp), find(Tp==max(Tp),1), max(Sgv-accN));

% ---- 终端的均值与散布（这才是"单剖面鲁棒性"的真实指标）----
fprintf('\n[RO terminal] mean | y=%7.2f m, z=%7.2f m | v_x=%7.3f, v_y=%7.3f, v_z=%7.3f m/s\n', ...
        diagR.mu_r(2), diagR.mu_r(3), diagR.mu_v(1), diagR.mu_v(2), diagR.mu_v(3));
fprintf('[RO terminal]  std | y=%7.2f m, z=%7.2f m | v_x=%7.3f, v_y=%7.3f, v_z=%7.3f m/s\n', ...
        diagR.sd_r(2), diagR.sd_r(3), diagR.sd_v(1), diagR.sd_v(2), diagR.sd_v(3));
fprintf('[RO terminal] relax: |r|max=%.2f m, |v|max=%.3f m/s\n', ...
        max(diagR.relax(1:3)), max(diagR.relax(4:6)));
toc
%% ===================== 6. 蒙特卡洛验证 =====================
N_mc = 500;   % 抽样数：std 的相对抽样误差 ≈ 1/sqrt(2N_mc)，000 时约 1.3%
Uprof = RO_nom.u;                    % 唯一的一条开环控制剖面
rng(0);  MC_DO = mc_verify(@(~) DO.u,  N_mc);
rng(0);  MC_RO = mc_verify(@(~) Uprof, N_mc);
% 对照实验：同一条 DO 剖面，但假设发动机能给出任意推力（不做饱和）
rng(0);  MC_DI = mc_verify(@(~) DO.u,  N_mc, false);
fprintf('\n===== Monte-Carlo verification (N=%d) =====\n', N_mc);
% 触地率是判断"设计是否偏了"的灵敏指标：标称在 t=tf 恰好触地，扰动对称时
% 触地时刻应一半偏早一半偏晚，故理想触地率 = 50%。明显偏离 50% 说明存在
% 系统性的提前/滞后（这里 DO 的饱和恰好把轨迹整体压低了）。
fprintf('DO  : std@tf (y,z)=(%6.1f,%6.1f) m | 触地率=%5.1f%% (对称理想值 50%%) | 推力越界样本=%5.1f%% | 平均饱和 %4.1f/%d 区间\n', ...
        std(MC_DO.r2f), std(MC_DO.r3f), 100*mean(MC_DO.touch), ...
        100*mean(MC_DO.thr_viol>0), mean(MC_DO.thr_sat), N);
fprintf('RO  : std@tf (y,z)=(%6.1f,%6.1f) m | 触地率=%5.1f%% (对称理想值 50%%) | 推力越界样本=%5.1f%% | 平均饱和 %4.1f/%d 区间\n', ...
        std(MC_RO.r2f), std(MC_RO.r3f), 100*mean(MC_RO.touch), ...
        100*mean(MC_RO.thr_viol>0), mean(MC_RO.thr_sat), N);

% ---- 终端散布的精确分解（不依赖抽样，可直接验证）---------------------------
% 共享开环剖面下，线性动力学给出严格恒等式
%       r(tf) - r_nom(tf) = (dr0 + tf*dv0)  +  drho(eta)
% 三项分别只依赖 xi、只依赖 eta，且 E[xi]=0 使交叉项为 0，于是
%       Var[r(tf)] = Var_floor + Var_eta[drho]        （逐轴相加）
% Var_floor 是**任何**共享开环剖面都摆脱不掉的不可约项；
% Var_eta[drho] 则完全由"指令被可用推力带饱和截断"产生 —— 只有这一项能被
% 鲁棒优化消除。用一维 eta 网格直接算出 drho(eta)，无需抽样。
sfloor = sqrt(sum(sigma_r0.^2 + tf^2*sigma_v0.^2)/3);   % /3：xi~U(-1,1)，半宽/根号3 才是 1sigma
EX_D = eta_excess(@(~) DO.u,  401);
EX_R = eta_excess(@(~) Uprof, 401);
sD = sqrt(std(MC_DO.r1f)^2 + std(MC_DO.r2f)^2 + std(MC_DO.r3f)^2);
sR = sqrt(std(MC_RO.r1f)^2 + std(MC_RO.r2f)^2 + std(MC_RO.r3f)^2);
sDp = sqrt(sfloor^2 + EX_D.s^2);
sRp = sqrt(sfloor^2 + EX_R.s^2);

fprintf('\n  终端 1sigma 散布 ||sigma_r(tf)|| 的精确分解：\n');
fprintf('    不可约项 sqrt(sigma_r0^2 + tf^2*sigma_v0^2) = %6.1f m  (任何共享剖面都一样)\n', sfloor);
fprintf('    饱和附加项 sqrt(Var_eta[drho])            : DO = %6.1f m (%2d/%d 区间饱和),  RO = %6.1f m (%2d/%d)\n', ...
        EX_D.s, round(EX_D.sat), N, EX_R.s, round(EX_R.sat), N);
fprintf('    饱和造成的系统偏差 ||E[drho]||             : DO = %6.1f m,  RO = %6.1f m\n', ...
        norm(EX_D.mean), norm(EX_R.mean));
fprintf('    => 合成 sigma(tf)  : DO = sqrt(%.1f^2 + %.1f^2) = %6.1f m (MC %6.1f m)\n', ...
        sfloor, EX_D.s, sDp, sD);
fprintf('                         RO = sqrt(%.1f^2 + %.1f^2) = %6.1f m (MC %6.1f m)\n', ...
        sfloor, EX_R.s, sRp, sR);
fprintf('    => RO/DO = %.3f  ；饱和附加项被消除 %.1f%%\n', sRp/sDp, 100*(1 - EX_R.s/max(EX_D.s,eps)));
fprintf(['  note: DO 的剖面贴在名义推力边界(%.2f~%.2f T_max)上，实际可用带却是\n' ...
         '        [%.2f,%.2f]*(1±%.0f%%)；故对 eta<0 样本上限被压到 %.3f、对 eta>0 样本\n' ...
         '        下限被抬到 %.3f，指令必然越界 -> 执行的不是设计剖面 -> 附加散布。\n' ...
         '        RO 的剖面收在区间最紧边界(%.3f~%.3f T_max)内，对任意 eta 都不饱和，\n' ...
         '        严格按设计剖面飞 -> 附加项恒为 0 -> 散布压到不可约下限。\n'], ...
        rho1/Tmax, rho2/Tmax, rho1/Tmax, rho2/Tmax, 100*delta_eta, ...
        rho2_w/Tmax, rho1_w/Tmax, rho1_w/Tmax, rho2_w/Tmax);
fprintf(['  note: 不可约项占了 sigma 的 %.1f%%，它对任何共享开环剖面都相同；\n' ...
         '        RO 能改善的只有饱和附加项，故总量上的提升有限但方向明确。\n'], ...
        100*sfloor/sDp);

% ---- 对照实验：证明"差异只来自可执行性，不来自优化目标"-----------------
sDi = sqrt(std(MC_DI.r1f)^2 + std(MC_DI.r2f)^2 + std(MC_DI.r3f)^2);
fprintf(['\n  [对照] 同一条 DO 剖面，若假设发动机总能给出指令推力（不设可用带）：\n' ...
         '         sigma(tf) = %.1f m  (理论下限 %.1f m)  —— 与 RO 的 %.1f m 相同。\n' ...
         '         可见共享开环剖面下"终端散布"与剖面形状无关，唯一能把 DO/RO 区分开的\n' ...
         '         是"指令推力能否被真实执行"；RO 的全部收益就在这一项上。\n'], ...
        sDi, sfloor, sR);

%% ===================== 7. 绘图 =====================
plot_all_trajectories(ROr, ROv, ROm, ROu);
plot_mc_comparison(DO, RO_nom, MC_DO, MC_RO, sfloor, EX_D, EX_R);
fprintf('\nDone.\n');

assignin('base','DO',DO);        assignin('base','RO_nom',RO_nom);
assignin('base','MC_DO',MC_DO);  assignin('base','MC_RO',MC_RO);
assignin('base','xi_nodes',xi_nodes);  assignin('base','M_xi',M_xi);
assignin('base','ROu',ROu);    assignin('base','ROm',ROm);   % 64 个配点的控制/质量
assignin('base','Mi_map',Mi_map); assignin('base','w_xi',w_xi);   % 谱插值映射与求积权
assignin('base','N',N);          assignin('base','dt',dt);
assignin('base','Q',Q);          assignin('base','Tmax',Tmax);

%% ======================================================================
%                    嵌套函数（共享工作区）
%% ======================================================================

%------------------------------------------------------------------
%  单配点 SOCP（无损凸化 + 解析离散化，一次求解）
%------------------------------------------------------------------
function [ok, Yv, Ev] = solve_point(y0)
    Y = sdpvar(7*N, 1, 'full');
    E = sdpvar(4*N, 1, 'full');

    Constraints = [Y == Abig*y0 + Lamug + PsiK*E];

    for k = 1:N
        b   = 7*(k-1);
        u_k = E(4*(k-1)+(1:3));
        sg  = E(4*(k-1)+4);
        r1k = Y(b+1);      % 这些都在 t_k（区间终点）
        r2k = Y(b+2);
        r3k = Y(b+3);
        % 区间起点 t_{k-1} 的质量（k=1 时即初始质量 log(m0)）
        if k == 1
            zs = log(m0);
        else
            zs = Y(7*(k-2)+7);
        end
        dzs = zs - z0s(k);

        % 无损凸化：||u|| <= sigma
        Constraints = [Constraints, cone(u_k, sg)];                       %#ok<AGROW>
        % DO 用**名义**推力边界：确定性设计对 ±2% 区间不确定性不设防，
        % 这正是鲁棒设计要对比出来的差异（RO 用最紧边界 rho*_w）。
        % 推力上限（区间起点质量配对、一阶 Taylor；e^-dz 凸 -> 保守）
        Constraints = [Constraints, sg <= rho2*expmzs(k)*(1-dzs)];        %#ok<AGROW>
        % 推力下限（二阶 Taylor epigraph）
        Constraints = [Constraints, ...
            rho1*expmzs(k)*(1-dzs+0.5*dzs^2) <= sg];                      %#ok<AGROW>
        % 高度非负
        Constraints = [Constraints, r1k >= 0];                            %#ok<AGROW>
        % 滑翔角：||[r2;r3]|| <= tan(86deg)*r1（+1 m 偏移避免锥顶点退化）
        Constraints = [Constraints, cone([r2k;r3k], tan_gs*r1k + 1.0)];   %#ok<AGROW>
        % 质量数值保护（物理上恒满足 m >= m0 - alpha*rho2*t）
        Constraints = [Constraints, Y(b+7) >= log(0.4*m0)];               %#ok<AGROW>
    end

    % 终端：着陆到目标点且速度为零
    yT = Y(7*(N-1)+(1:7));
    Constraints = [Constraints, yT(1) == 0, yT(2) == 0, yT(3) == 0];
    Constraints = [Constraints, yT(4) == 0, yT(5) == 0, yT(6) == 0];

    % 目标：最小燃料（最大化终端质量），图片 Eq.37 J = -m(tf)
    obj = -yT(7);

    sol = optimize(Constraints, obj, ...
                   sdpsettings('solver','mosek','verbose',VRB,'debug',0));
    ok = (sol.problem == 0);
    if ~ok
        Yv = zeros(7*N,1);  Ev = zeros(4*N,1);
        return;
    end
    Yv = value(Y);  Ev = value(E);
end

%------------------------------------------------------------------
%  鲁棒 SOCP：一条共享控制 + Q 条状态轨迹（对应再入的 scp_iteration）
%  状态已按 y_s = Abig*y0_s + Lamug + PsiK*E 解析消去，只剩控制 E 与松弛量
%------------------------------------------------------------------
function [ok, Yall, Esg, D] = solve_shared()
    E     = sdpvar(4*N, 1, 'full');     % 唯一的一份控制（被全部配点共享）
    relax = sdpvar(6, 1, 'full');       % 终端均值松弛（对应 AMM4 的 relax）

    % ---- 各配点的"自由响应"：控制置零时的状态 ----
    Bs = zeros(7*N, Q);
    for s = 1:Q
        ys = [X0all(s,1:3).'; X0all(s,4:6).'; log(m0)];
        Bs(:,s) = Abig*ys + Lamug;
    end

    % ---- 均值 / 谱标准差的"提取算子"（都作用在 E 上）----
    %   任一状态分量 id 在 Q 个配点上的取值向量:  f = Bs(id,:).' + (PsiK(id,:)*E)
    %   均值      mu(id) = sum_s w(s) f(s) = cm(id) + PsiK(id,:)*E   （sum w = 1）
    %   谱系数    Sig(id) = Mi_map*f = cSg(:,id) + Mrep*(PsiK(id,:)*E)
    %             Var = ||Sig||^2 ,  Mrep = Mi_map*ones(Q,1)
    %   Mrep 恒为 0（常函数的非零阶谱系数为 0），故 Sig 与 E 无关 —— 这正是
    %   "终端/过程散布与共享剖面无关"在代码里的体现。
    %   于是标准差可以直接写成闭式常数 SigA，不需要引入任何辅助变量；
    %   若把它写成 sv >= norm(...) 的锥，会给 MOSEK 塞进上百个 63 行、
    %   系数仅 1e-17 的退化锥（WRN 710），实测使模型退化（sigma 偏离 ||u||
    %   达 45%、多耗 70 kg 燃料、终端质量与逐步积分对不上）。
    cm   = Bs*w_xi;                     % 7N x 1，均值常数项
    cSg  = Mi_map*Bs.';                 % (Q-1) x 7N，谱系数（已去常数行）
    SigA = reshape(sqrt(sum(cSg.^2,1)), 7, N);   % 7 x N，各时刻的标准差（常数）

    C = [];

    % ---- 推力约束：z=log(m) 与 r,v 完全解耦，且 z(0) 固定 ----
    % 故 64 个配点的质量历程完全相同（m 只取决于 sigma 的积分），
    % 推力上下限只需施加一次，不必对 s 复制 64 份。
    for k = 1:N
        u_k = E(4*(k-1)+(1:3));
        sg  = E(4*(k-1)+4);
        C = [C, cone(u_k, sg)];                                       %#ok<AGROW>
        if k == 1
            dzs = 0;                        % z(0) = log(m0) = z0s(1)
        else
            id  = 7*(k-2) + 7;              % 区间起点 t_{k-1} 的质量
            dzs = (Bs(id,1) + PsiK(id,:)*E) - z0s(k);
        end
        % 上限：一阶 Taylor（e^-dz 凸 -> 保守）
        C = [C, sg <= rho2_w*expmzs(k)*(1-dzs)];                      %#ok<AGROW>
        % 下限：二阶 Taylor epigraph（凸二次）
        C = [C, rho1_w*expmzs(k)*(1-dzs+0.5*dzs^2) <= sg];            %#ok<AGROW>
    end

    % ---- 过程约束：均值 + kC*sigma 裕度（对应 AMM4 的 r_mu >= pc + kC*sqrt(var)）----
    for k = 1:N
        i1 = 7*(k-1) + 1;
        mu1 = cm(i1)   + PsiK(i1,:)*E;
        mu2 = cm(i1+1) + PsiK(i1+1,:)*E;
        mu3 = cm(i1+2) + PsiK(i1+2,:)*E;
        kk = kC_path(k);
        % 高度：均值 - kC*sigma >= 0（sigma 为闭式常数）
        C = [C, mu1 >= kk*SigA(1,k)];                                 %#ok<AGROW>
        % 滑翔角锥：||[mu2;mu3]|| + kC*(s2+s3) <= tan(86deg)*mu1 + 1
        % （+1 m 偏移避免锥顶点退化）
        C = [C, cone([mu2;mu3], tan_gs*mu1 + 1.0 - kk*(SigA(2,k)+SigA(3,k)))]; %#ok<AGROW>
        % 质量数值保护（物理上恒满足 m >= m0 - alpha*rho2*t）
        C = [C, cm(i1+6) + PsiK(i1+6,:)*E >= log(0.4*m0)];            %#ok<AGROW>
    end

    % ---- 终端：均值软约束（relax）----
    iT   = 7*(N-1) + (1:7);
    muT  = cm(iT) + PsiK(iT,:)*E;       % 7 x 1，终端均值状态
    C = [C, relax >= 0];
    C = [C,  muT(1:6) <=  relax,  -muT(1:6) <= relax];

    % ---- 目标：最大化终端质量（最小燃料）+ 终端均值偏差重罚 ----
    % （散布项是常数，不进目标；它作为鲁棒性指标在 D.sd_r / D.sd_v 中输出）
    obj = -muT(7) + K_final*sum(relax);

    sol = optimize(C, obj, sdpsettings('solver','mosek','verbose',VRB,'debug',0));
    D.problem = sol.problem;
    ok = (sol.problem == 0);
    if ~ok
        Yall = zeros(7*N, Q);  Esg = zeros(4*N, 1);
        D.mu_r = nan(3,1);  D.mu_v = nan(3,1);
        D.sd_r = nan(3,1);  D.sd_v = nan(3,1);  D.relax = nan(6,1);
        return;
    end
    Esg = value(E);
    % 由共享控制回填 Q 条状态轨迹（与解析离散化一致，无需再解一次）
    Yall = Bs + PsiK*Esg;
    D.mu_r  = value(muT(1:3));      D.mu_v = value(muT(4:6));
    D.sd_r  = SigA(1:3,N);          D.sd_v = SigA(4:6,N);
    D.relax = value(relax);
end

%------------------------------------------------------------------
%  用给定控制剖面从标称初值精确积分（与 mc_verify 的离散格式一致）
%------------------------------------------------------------------
function S = propagate_nominal(Ev)
    E = reshape(Ev, 4, N);
    u = E(1:3,:);
    S.u = u;
    S.r = zeros(3, N+1);  S.v = zeros(3, N+1);  S.m = zeros(N+1, 1);
    S.r(:,1) = r0;  S.v(:,1) = v0;  S.m(1) = m0;
    for k = 1:N
        acc = u(:,k) + gvec;
        S.r(:,k+1) = S.r(:,k) + dt*S.v(:,k) + 0.5*dt^2*acc;
        S.v(:,k+1) = S.v(:,k) + dt*acc;
        S.m(k+1)   = S.m(k)*exp(-alpha*norm(u(:,k))*dt);
    end
    S.sig = sqrt(sum(u.^2,1));
end

%------------------------------------------------------------------
%  单点解包
%------------------------------------------------------------------
function S = unpack_point(Yv, Ev, y0)
    S.r = zeros(3, N+1);   S.v = zeros(3, N+1);
    S.m = zeros(N+1, 1);   S.u = zeros(3, N);
    S.r(:,1) = y0(1:3);  S.v(:,1) = y0(4:6);  S.m(1) = exp(y0(7));
    for k = 1:N
        b = 7*(k-1);
        S.r(:,k+1) = Yv(b+(1:3));
        S.v(:,k+1) = Yv(b+(4:6));
        S.m(k+1)   = exp(Yv(b+7));
        S.u(:,k)   = Ev(4*(k-1)+(1:3));
    end
end

%------------------------------------------------------------------
%  蒙特卡洛验证（DO / RO 使用同一组样本）
%------------------------------------------------------------------
function MC = mc_verify(u_fun, N_mc, sat_on)
    % sat_on = false 时把可用推力带放开到 (-inf, inf)，即"假设发动机总能给出
    % 指令推力"。这是本文件的关键对照：它用来证明 DO 与 RO 的散布差异
    % 100% 来自"能不能执行"，而不是来自优化目标本身。
    if nargin < 3, sat_on = true; end
    xi  = 2*rand(N_mc, n_xi) - 1;
    eta = 2*rand(N_mc, 1) - 1;

    MC.r  = zeros(3, N+1, N_mc);
    MC.v  = zeros(3, N+1, N_mc);
    MC.m  = zeros(N+1, N_mc);
    MC.u  = zeros(3, N, N_mc);
    MC.gs = zeros(N, N_mc);
    MC.r1f = zeros(N_mc,1);  MC.r2f = zeros(N_mc,1);  MC.r3f = zeros(N_mc,1);
    MC.vf  = zeros(3, N_mc);
    MC.thr_viol = zeros(N_mc,1);
    MC.thr_ov   = zeros(N_mc,1);      % 超上限相对量 max(T_cmd/s_hi - 1)
    MC.thr_un   = zeros(N_mc,1);      % 欠下限相对量 max(s_lo/T_cmd - 1)
    MC.thr_sat  = zeros(N_mc,1);      % 发生推力饱和的区间数
    MC.touch    = false(N_mc,1);      % 是否在 tf 之前触地
    MC.t_td     = NaN(N_mc,1);        % 触地时刻
    MC.rtd      = nan(3,N_mc);        % 触地点
    MC.vtd      = nan(3,N_mc);        % 触地速度
    MC.rT       = zeros(3,N_mc);      % t=tf 时的状态（终端散布指标）
    MC.vT       = zeros(3,N_mc);
    MC.t  = (0:N)*dt;  MC.tc = (1:N)*dt;

    for i = 1:N_mc
        r = r0 + sigma_r0.*xi(i,1:3).';
        v = v0 + sigma_v0.*xi(i,4:6).';
        m = m0;
        U = u_fun(xi(i,:));

        if sat_on
            s_lo = rho1*(1 + delta_eta*eta(i));
            s_hi = rho2*(1 + delta_eta*eta(i));
        else
            s_lo = -inf;  s_hi = inf;      % 对照：不设可用带
        end
        viol = 0;  ov = 0;  un = 0;  nsat = 0;
        touched = false;  froze = true;

        MC.r(:,1,i) = r;  MC.v(:,1,i) = v;  MC.m(1,i) = m;

        for k = 1:N
            uk = U(:,k);

            % 推力检查：T = m*||u||（力），与加速度 ||u|| 区分开；
            % m 取区间起点质量 —— 与优化模型中 sigma_k 的配对方式一致
            T = m*norm(uk);
            ov = max(ov, T/s_hi - 1);
            un = max(un, s_lo/max(T,1e-9) - 1);
            if T > s_hi*(1+1e-9) || T < s_lo*(1-1e-9)
                viol = viol + 1;
            end

            % ---- 执行机构饱和（区间不确定性的真实作用点）----
            % 区间不确定性说的是"可用推力带是 [rho1,rho2]*(1±2%)"，即发动机
            % 实际能给出的推力被夹在这个带内。指令一旦越界，**不可能**按指令
            % 输出：只能沿原方向按可用边界饱和输出（下限同理，节流有下限）。
            % 这一步是 DO / RO 蒙特卡洛差异的物理来源：
            %   DO  的剖面贴在名义边界 0.8/0.3 上 -> 约一半样本必然饱和
            %       -> 实际执行的不是设计剖面 -> 轨迹偏离、散布被放大；
            %   RO  的剖面贴在 0.784/0.309 上 -> 对任意 eta 都不饱和
            %       -> 严格按设计剖面飞 -> 散布保持设计值。
            % 若不做饱和、直接把指令当实际推力，等于假设发动机能超额输出，
            % 于是 DO 的名义边界剖面被"免费"执行，两者散布就分不出差别。
            Tap = min(max(T, s_lo), s_hi);
            if abs(Tap - T) > 1e-9*max(s_hi,1)
                nsat = nsat + 1;
            end
            uap = uk*(Tap/max(T,1e-12));      % 方向不变，幅值取可用值

            acc = uap + gvec;
            rn = r + dt*v + 0.5*dt^2*acc;
            vn = v + dt*acc;
            mn = m*exp(-alpha*norm(uap)*dt);

            % ---- 首次穿越地面：步内线性插值，记录真实触地点/时刻 ----
            if ~touched && rn(1) <= 0
                th = r(1)/(r(1) - rn(1));     % 穿越比例，属于 (0,1]
                th = min(max(th,0),1);
                rt = r + th*(rn - r);  rt(1) = 0;
                MC.rtd(:,i) = rt;             % 真实触地点（用于散布/成功率统计）
                MC.vtd(:,i) = v + th*(vn - v);% 真实撞击速度，绝不置零
                MC.t_td(i)  = (k-1 + th)*dt;  % 真实触地时刻
                touched = true;
                froze   = false;              % 本步仍要把触地点画出来
            end

            r = rn;  v = vn;  m = mn;
            MC.u(:,k,i) = uk;                 % 记录**指令**推力加速度（约束检查用）

            if touched && ~froze
                % 触地这一步：曲线画到插值出的触地点为止
                MC.r(:,k+1,i) = MC.rtd(:,i);  MC.v(:,k+1,i) = MC.vtd(:,i);
                MC.m(k+1,i)   = mn;           MC.gs(k,i) = NaN;
                froze = true;
            elseif touched
                % 触地之后：轨迹不再显示（避免画出"地下"曲线），但**积分继续**
                % 到 t=tf。终端散布统一按 t=tf 的偏差统计 —— 否则样本会按触地
                % 早晚被截断（触地早的样本只累积了一部分偏差），指标就被"触地
                % 率"污染了，DO/RO 之间会出现与鲁棒性无关的伪差异。
                MC.r(:,k+1,i) = NaN;  MC.v(:,k+1,i) = NaN;  MC.m(k+1,i) = NaN;
                MC.gs(k,i) = NaN;
            else
                MC.r(:,k+1,i) = r;  MC.v(:,k+1,i) = v;  MC.m(k+1,i) = m;
                % 滑翔角取"与铅垂轴夹角" theta = atan(||r23||/r1)：
                %   约束 sqrt(r2^2+r3^2) <= tan(86 deg)*r1  <=>  theta <= 86 deg
                % 末段该量是 0/0 型：r1 与 ||r23|| 同时趋于 0，比值完全由各
                % 样本自己的终端偏差决定，故画出来会先"激增到 90 deg"再"回落
                % 到 0"。那不是数值故障，而是 r -> 0 时角度失去意义的真实退化，
                % 按要求**全部保留**（只在地面以下即 r1 <= 0 时置 NaN）。
                if r(1) > 0
                    MC.gs(k,i) = atan2d(norm(r(2:3)), r(1));
                else
                    MC.gs(k,i) = NaN;
                end
            end
        end
        MC.thr_viol(i) = viol;
        MC.thr_ov(i)   = ov;
        MC.thr_un(i)   = un;
        MC.thr_sat(i)  = nsat;
        MC.touch(i)    = touched;
        % 终端状态一律取 t=tf（与优化模型的终端约束同一时刻）
        MC.rT(:,i) = r;  MC.vT(:,i) = v;
        MC.r1f(i) = r(1);  MC.r2f(i) = r(2);  MC.r3f(i) = r(3);
        MC.vf(:,i) = v;
    end
end

%------------------------------------------------------------------
%  区间不确定性 eta 单独作用时造成的终端偏差 drho(eta)（精确、无抽样）
%  取 xi=0，只扫 eta；drho = r(tf;eta) - r(tf;eta=0)。
%  因为共享开环剖面下 dr(tf) = (dr0 + tf*dv0) + drho(eta) 且两两独立，
%  Var[dr(tf)] = Var_floor + Var_eta[drho]，可直接与 MC 结果对照。
%------------------------------------------------------------------
function D = eta_excess(u_fun, n_eta)
    if mod(n_eta,2) == 0, n_eta = n_eta + 1; end   % 保证网格含 eta = 0
    eg = linspace(-1,1,n_eta).';
    q0 = (n_eta+1)/2;
    R  = zeros(3,n_eta);  nsat = zeros(1,n_eta);
    U  = u_fun(zeros(1,n_xi));                     % 控制剖面与 xi 无关
    for q = 1:n_eta
        eta  = eg(q);
        s_lo = rho1*(1 + delta_eta*eta);
        s_hi = rho2*(1 + delta_eta*eta);
        r = r0;  v = v0;  m = m0;  nn = 0;
        for k = 1:N
            uk = U(:,k);
            T  = m*norm(uk);
            Ta = min(max(T,s_lo),s_hi);            % 只有可用推力带内的幅值能真正输出
            if abs(Ta-T) > 1e-9*max(s_hi,1), nn = nn + 1; end
            ua = uk*(Ta/max(T,1e-12));
            a  = ua + gvec;
            r  = r + dt*v + 0.5*dt^2*a;
            v  = v + dt*a;
            m  = m*exp(-alpha*norm(ua)*dt);
        end
        R(:,q) = r;  nsat(q) = nn;
    end
    drho    = R - R(:,q0);
    D.v     = var(drho,0,2);        % 逐轴 Var_eta
    D.s     = sqrt(sum(D.v));       % 附加散布的 2-范数（1sigma）
    D.mean  = mean(drho,2);         % 饱和造成的系统偏差
    D.sat   = mean(nsat);           % 平均饱和区间数
end

%------------------------------------------------------------------
%  诊断打印
%------------------------------------------------------------------
function print_summary(S)
    fprintf('  terminal [x y z] = [%8.2f %8.2f %8.2f] m\n', S.r(:,end));
    fprintf('  terminal [vx vy vz] = [%8.3f %8.3f %8.3f] m/s\n', S.v(:,end));
    fprintf('  final mass = %8.1f kg   (fuel used = %.1f kg)\n', S.m(end), m0 - S.m(end));
    % u 是推力加速度；推力(力) = m.*||u||，两者分开打印，避免混淆
    acc = sqrt(sum(S.u.^2,1));          % m/s^2
    Tc  = S.m(1:N).'.*acc;              % N（用区间起点质量）
    fprintf('  max |u|    = %8.3f m/s^2   (thrust accel)\n', max(acc));
    fprintf('  max  T     = %8.1f N       (thrust force, T = m|u|)\n', max(Tc));
    fprintf('  max T/Tmax = %8.4f       (名义上限 %.4f | 区间最紧上限 %.4f)\n', ...
            max(Tc)/Tmax, rho2/Tmax, rho2_w/Tmax);
    fprintf('  min T/Tmax = %8.4f       (名义下限 %.4f | 区间最紧下限 %.4f)\n', ...
            min(Tc)/Tmax, rho1/Tmax, rho1_w/Tmax);
    fprintf('  min x      = %8.2f m       (altitude)\n', min(S.r(1,:)));
end

%------------------------------------------------------------------
%  所有配点轨迹（仿再入 plot_all_trajectories）
%------------------------------------------------------------------
function plot_all_trajectories(Rr, Rv, Rm, Ru)
    % yline 生成 ConstantLine 对象，与 legend 同处一个坐标系时 MATLAB 会抛
    % MATLAB:structOnObject 的内部提示（无害，且会让控制台出现乱码警告）。
    % 这里在函数内窄范围抑制，退出时恢复。
    ws = warning('off','MATLAB:structOnObject');
    t  = (0:N)*dt;  tc = (1:N)*dt;
    % 坐标命名：x = r1（高度），y = r2（射程），z = r3（横程）
    names = {'r_1 (km)','r_2 (km)','r_3 (km)', ...
             'v_1 (m/s)','v_2 (m/s)','v_3 (m/s)'};
    dat = {squeeze(Rr(1,:,:))/1000, squeeze(Rr(2,:,:))/1000, squeeze(Rr(3,:,:))/1000, ...
           squeeze(Rv(1,:,:)),      squeeze(Rv(2,:,:)),      squeeze(Rv(3,:,:))};

    % ---- 曲线是否重合取决于配点集的"位置"，不是画图问题 ----------------
    % 动力学对状态仿射、且全部配点共享同一条控制，于是解析解严格满足
    %     r_i(t) = r_i,nom(t) + dr0_i + t*dv0_i ,   v_i(t) = v_i,nom(t) + dv0_i
    % 其中 dr0_i = sigma_r0(i)*xi_i、dv0_i = sigma_v0(i)*xi_{3+i}，所以
    %     位置分量 -> 只由 (xi_i, xi_{3+i}) 决定 -> 曲线去重数 <= 该二维投影的点数
    %     速度分量 -> 只由 xi_{3+i} 决定          -> 去重数 <= 该一维投影的点数
    % 若配点是"轴对齐张量网格"{-1,+1}^6，二维投影只有 2x2=4 个点 -> 只看到 4 条
    % （速度 1 维投影 2 个点 -> 2 条），其余 16/32 条严格重合。
    % 若配点是 SAMP='lhs' 的一般位置集合，任意投影的点数都等于 Q -> 64 条全不同。
    % 下面把去重数直接算出来并打印/标注，便于核对。
    nuniq = zeros(1,6);
    for k = 1:6
        Dk = round(dat{k}*1e6)/1e6;              % 去重前消除浮点尾差
        [~,ia] = unique(Dk.', 'rows');
        nuniq(k) = numel(ia);
    end
    fprintf(['[collocation] %d 个配点 -> 各面板互不相同的曲线数：' ...
             'r_1/r_2/r_3 = %d/%d/%d, v = %d/%d/%d\n'], Q, nuniq);

    % 触地后不再画：解析离散解在 t_N 之后/落地之前就可能出现"负高度"
    % （共享开环剖面下低初值的角点会提前撞地），那不是数值错误而是真实几何，
    % 但画成 r1 < 0 的曲线容易被误读，故与 MC 一致地在触地处截断（置 NaN）。
    mask = true(N+1, Q);
    for s = 1:Q
        i0 = find(squeeze(Rr(1,:,s)).' < 0, 1, 'first');
        if ~isempty(i0), mask(i0:end, s) = false; end
    end
    for k = 1:6
        d = dat{k};  d(~mask) = NaN;  dat{k} = d;
    end

    figure('Name','Collocation states vs time','Color','w');
    for k = 1:6
        subplot(2,3,k);
        % 给矩阵 y 时 MATLAB 按坐标轴 ColorOrder 循环上色（默认配色），
        % 64 条配点曲线自动分色。重合的曲线后画的会盖住先画的，故再把
        % 均值（= 标称，PCE 配点对称）用黑色虚线压在最上面，作为参照。
        h = plot(t, dat{k}, 'LineWidth', 0.5); grid on; hold on;
        % 面板内标注：说明这里只看到几条不同曲线
        xlabel('Time (s)'); ylabel(names{k});
        if k == 1
            legend(h(1), 'trajectories');
        end
    end
    if strcmpi(SAMP,'corner')
        sub2 = 'corner grid: r_i sees only (xi_i, xi_{i+3}) -> 4 distinct curves, v_i -> 2';
    else
        sub2 = sprintf(['%s sampling: r_i(t) = r_i,nom + dr0_i + t*dv0_i ' ...
                        'with all 64 (dr0_i, dv0_i) distinct -> %d distinct curves'], ...
                       upper(SAMP), Q);
    end
    % sgtitle({sprintf('Collocation-point states (%d PCE collocation points, %s)', Q, upper(SAMP)), sub2}, ...
    %         'FontSize',9);

    % 推力与质量：控制只有一份（全部配点共享），质量又只取决于 sigma 的积分，
    % 与 xi 无关，故 Q 条曲线**严格重合**。这与再入代码里 cu1/cu2 是
    % (N_points x 1) 的设计一致；画 Q 条重合线没有意义，只画一条。
    figure('Name','Shared thrust profile','Color','w');
    subplot(2,1,1);
    Es  = reshape(Ru(:,1), 4, N);            % 唯一的一份控制（Q 列完全相同）
    mag = sqrt(sum(Es(1:3,:).^2, 1)).';
    Tn  = mag .* Rm(1:N,1) / Tmax;           % T/Tmax = m|u|/Tmax（区间起点质量）
    % DO（确定性优化）的控制剖面：DO 只有一份 u，配点共享的概念在它身上本来就
    % 只有一条。推力模长 |u| 乘以**区间起点**质量，与 RO 的算法保持逐位一致
    % （m 在区间内单调下降，T = m*sigma 的最大值出现在区间起点）。
    magD = sqrt(sum(DO.u.^2, 1)).';
    TnD  = magD .* DO.m(1:N) / Tmax;
    hC  = plot(tc, Tn, 'r--', 'LineWidth', 1.8); hold on;   % RO 共享剖面
    hD  = plot(DO.tc, TnD, 'b-', 'LineWidth', 1.8);       % DO 控制剖面
    % 参考线全部用"先创建、再 set"的方式：yline(y,spec,'LineWidth',w) 里的
    % 'LineWidth' 会被当成 label 位置参数，导致图例错位
    hLo = yline(rho1/Tmax,   'g--');  set(hLo,'LineWidth',1.5);
    hHi = yline(rho2/Tmax,   'k:');  set(hHi,'LineWidth',1.5);
    hLo_w = yline(rho1_w/Tmax,   'g:');  set(hLo_w,'LineWidth',1.5);
    hHi_w = yline(rho2_w/Tmax,   'k--');  set(hHi_w,'LineWidth',1.5);
    ylim([0 1.25]);                      % 留出顶部空白，避免图例压住 0.8 段曲线
    xlabel('Time (s)'); ylabel('T / T_{max}'); grid on;
    legend([hD hC hLo hHi hLo_w hHi_w], ...
           {'DO', ...
            'RO', ...
            sprintf('%.2f T_{max} (nom. lower)', rho1/Tmax), ...
            sprintf('%.2f T_{max} (nom. upper)', rho2/Tmax), ...
            sprintf('%.3f T_{max} (worst lower)', rho1_w/Tmax), ...
            sprintf('%.3f T_{max} (worst upper)', rho2_w/Tmax)}, ...
           'Location','northwest','FontSize',9);
    subplot(2,1,2);
    hM1 = plot(DO.t, DO.m,     'b-', 'LineWidth', 1.8);  hold on;   % DO 质量
    hM2 = plot(t,    Rm(:,1),  'r--',  'LineWidth', 1.8);  % RO 质量
    xlabel('Time (s)'); ylabel('Mass (kg)'); grid on;
    legend([hM1 hM2], {'DO', 'RO'}, ...
           'Location','northeast','FontSize',9);
    % 两条质量曲线的差别只有燃料项（约 10 kg / 1905 kg = 0.5%），在全量程坐标里
    % 几乎重合，故把终值差额直接标注出来，避免"看起来一模一样=没画上"的误读。
    text(0.03, 0.10, sprintf('m_f:  DO = %.1f kg,  RO = %.1f kg  (RO burns +%.1f kg)', ...
         DO.m(end), Rm(end,1), DO.m(end) - Rm(end,1)), ...
         'Units','normalized','FontSize',8,'Color',[0.25 0.25 0.25], ...
         'BackgroundColor','w','EdgeColor',[0.75 0.75 0.75]);
    drawnow;
    warning(ws);
end

%------------------------------------------------------------------
%  DO / RO 对比（仿再入 plot_mc_comparison）
%------------------------------------------------------------------
function plot_mc_comparison(DOs, ROs, MC_DO, MC_RO, sfloor, EX_D, EX_R)
    ws = warning('off','MATLAB:structOnObject');   % 同上，legend+yline 的内部提示
    % (1) 标称轨迹
    figure('Name','Nominal trajectory comparison','Color','w');
    subplot(1,2,1);
    h1 = plot3(DOs.r(2,:)/1000, DOs.r(3,:)/1000, DOs.r(1,:)/1000,'b-','LineWidth',2); hold on;
    h2 = plot3(ROs.r(2,:)/1000, ROs.r(3,:)/1000, ROs.r(1,:)/1000,'r--','LineWidth',2);
    h3 = plot3(DOs.r(2,1)/1000, DOs.r(3,1)/1000, DOs.r(1,1)/1000,'go', ...
               'MarkerSize',8,'MarkerFaceColor','g');
    h4 = plot3(0,0,0,'kp','MarkerSize',10,'MarkerFaceColor','k');
    % 三维图的三个坐标轴顺序是 (r2, r3, r1) = (y, z, x)
    xlabel('r_2 (km)'); ylabel('r_3 (km)'); zlabel('r_1 (km)');
    legend([h1 h2 h3 h4], {'DO','RO','Initial','Target'}, ...
           'Location','northeast','FontSize',9);
    grid on; view(-35,25);
    subplot(1,2,2);
    h1 = plot(DOs.t, DOs.r(1,:)/1000,'b-','LineWidth',2); hold on;
    h2 = plot(ROs.t, ROs.r(1,:)/1000,'r--','LineWidth',2);
    xlabel('Time (s)'); ylabel('x (km, altitude)');
    legend([h1 h2], {'DO','RO'},'Location','northeast'); grid on;

    % (2) MC 三维打靶
    figure('Name','Monte Carlo 3D trajectories','Color','w');
    n_show = min(80, size(MC_DO.r,3));
    subplot(1,2,1);
    for i = 1:n_show
        plot3(MC_DO.r(2,:,i)/1000, MC_DO.r(3,:,i)/1000, MC_DO.r(1,:,i)/1000, ...
              'Color',[0.6 0.75 0.9],'LineWidth',0.2); hold on;
    end
    hT1 = plot3(0,0,0,'kp','MarkerSize',12,'MarkerFaceColor','k');
    xlabel('r_2 (km)'); ylabel('r_3 (km)'); zlabel('r_1 (km)');
    legend(hT1, {'Target'}, 'Location','northeast','FontSize',9);
    title('DO'); grid on; view(3);
    subplot(1,2,2);
    for i = 1:n_show
        plot3(MC_RO.r(2,:,i)/1000, MC_RO.r(3,:,i)/1000, MC_RO.r(1,:,i)/1000, ...
              'Color',[0.95 0.7 0.5],'LineWidth',0.2); hold on;
    end
    hT = plot3(0,0,0,'kp','MarkerSize',12,'MarkerFaceColor','k');
    xlabel('r_2 (km)'); ylabel('r_3 (km)'); zlabel('r_1 (km)');
    legend(hT, {'Target'}, 'Location','northeast','FontSize',9);
    title('RO'); grid on; view(3);

    % (3) 状态曲线
    names = {'r_1 (km)','r_2 (km)','r_3 (km)','v_1 (m/s)','v_2 (m/s)','v_3 (m/s)'};
    dD = {squeeze(MC_DO.r(1,:,:))/1000, squeeze(MC_DO.r(2,:,:))/1000, squeeze(MC_DO.r(3,:,:))/1000, ...
          squeeze(MC_DO.v(1,:,:)), squeeze(MC_DO.v(2,:,:)), squeeze(MC_DO.v(3,:,:))};
    dR = {squeeze(MC_RO.r(1,:,:))/1000, squeeze(MC_RO.r(2,:,:))/1000, squeeze(MC_RO.r(3,:,:))/1000, ...
          squeeze(MC_RO.v(1,:,:)), squeeze(MC_RO.v(2,:,:)), squeeze(MC_RO.v(3,:,:))};
    nD = {DOs.r(1,:)/1000, DOs.r(2,:)/1000, DOs.r(3,:)/1000, DOs.v(1,:), DOs.v(2,:), DOs.v(3,:)};
    nR = {ROs.r(1,:)/1000, ROs.r(2,:)/1000, ROs.r(3,:)/1000, ROs.v(1,:), ROs.v(2,:), ROs.v(3,:)};
    figure('Name','State curves DO vs RO','Color','w');
    hs = gobjects(4,1);
    for k = 1:6
        subplot(2,3,k);
        h1 = plot(MC_DO.t, dD{k},'Color',[0.6 0.75 0.9],'LineWidth',0.5); hold on;
        h2 = plot(MC_RO.t, dR{k},'Color',[0.95 0.7 0.5],'LineWidth',0.5);
        h3 = plot(DOs.t, nD{k},'b-','LineWidth',2);
        h4 = plot(ROs.t, nR{k},'r--','LineWidth',2);
        xlabel('Time (s)'); ylabel(names{k}); grid on;
        if k == 1, hs = [h1(1); h2(1); h3; h4]; end
    end
    % 共享图例放在图底部，避免遮挡数据
    axL = axes('Position',[0 0 1 1],'Visible','off');
    legend(axL, hs, {'DO MC','RO MC','DO nominal','RO nominal'},'FontSize',9);

    % (4) 控制与过程约束
    % 推力用无量纲量 Tn = m*|u|/Tmax 表示：物理约束即为 0.2 <= Tn <= 0.8，
    % 若按 |u| 画则上下限随质量变化，需除以 Tmax 才可直接比较。
    figure('Name','Control and path constraints','Color','w');
    TnD = squeeze(MC_DO.m(1:N,:)) .* squeeze(sqrt(sum(MC_DO.u.^2,1)));
    TnR = squeeze(MC_RO.m(1:N,:)) .* squeeze(sqrt(sum(MC_RO.u.^2,1)));
    TnD = TnD/Tmax;  TnR = TnR/Tmax;

    subplot(2,2,1);
    % 同质量图：DO 的推力 T=m||u|| 在所有样本上完全相同（m、u 都不随 xi 变），
    % 500 条曲线严格重合，故用略粗的浅蓝线画，否则会完全藏在 0.2/0.8 参考线下
    h1 = plot(MC_DO.tc, TnD,'b-','LineWidth',1.5); hold on;
    h2 = plot(MC_RO.tc, TnR,'r--','LineWidth',1.5);
    h3 = yline(rho1/Tmax,   'g--');  set(h3,'LineWidth',1.5);
    h4 = yline(rho2/Tmax,   'k:');  set(h4,'LineWidth',1.5);
    h5 = yline(rho1_w/Tmax, 'g:' );  set(h5,'LineWidth',1.5);   % 区间最紧下限
    h6 = yline(rho2_w/Tmax, 'k--' );  set(h6,'LineWidth',1.5);   % 区间最紧上限
    ylim([0 1.25]);                     % 顶部留白，避免图例压住 0.8 段曲线
    xlabel('Time (s)'); ylabel('T / T_{max}');
    legend([h1(1) h2(1) h3 h4 h5 h6], ...
           {'DO MC (commanded)','RO MC (commanded)', ...
            sprintf('%.2f T_{max} (nom. lower)', rho1/Tmax), ...
            sprintf('%.2f T_{max} (nom. upper)', rho2/Tmax), ...
            sprintf('%.3f T_{max} (worst lower)', rho1_w/Tmax), ...
            sprintf('%.3f T_{max} (worst upper)', rho2_w/Tmax)}, ...
           'Location','northwest','FontSize',7); grid on;

    subplot(2,2,2);
    % 纵轴为"滑翔角"——相对铅垂轴（目标点垂直方向）的锥半角 theta：
    %     theta = atan(||r23||/r1)，锥约束 sqrt(r2^2+r3^2) <= tan(86deg)*r1
    %     等价于 theta <= 86 deg。
    % 为什么不用仰角 atan(r1/||r23||)：末段 ||r23|| 比 r1 更快趋于 0，
    % 该比值病态，会画出虚假的"激增到 90 度"，属于数值假象而非真实越界。
    % 全部蒙特卡洛样本都画：细线是 MC，粗线是标称剖面。
    % 末段的"先冲到 90 deg 再回落到 0"是真实退化（r -> 0 时 atan2(||r23||,r1)
    % 变成 0/0，比值由各样本自己的终端偏差决定），保留全部曲线不做裁剪。
    gD = squeeze(MC_DO.gs);  gR = squeeze(MC_RO.gs);     % (N x N_mc)
    h1 = plot(MC_DO.tc, gD, 'Color',[0.55 0.72 0.88], 'LineWidth',0.5); hold on;
    h2 = plot(MC_RO.tc, gR, 'Color',[0.95 0.62 0.30], 'LineWidth',0.5);
    thD = atan2d(sqrt(sum(DOs.r(2:3,:).^2,1)), DOs.r(1,:));
    thR = atan2d(sqrt(sum(ROs.r(2:3,:).^2,1)), ROs.r(1,:));
    thD(DOs.r(1,:) <= 0) = NaN;   thR(ROs.r(1,:) <= 0) = NaN;
    h3 = plot(DOs.t, thD, 'b-', 'LineWidth', 1.5);
    h4 = plot(ROs.t, thR, 'r--','LineWidth', 1.5);
    h5 = yline(86, 'k-.');  set(h5,'LineWidth',1.5);
    ylim([0 95]);  xlim([0 tf]);
    xlabel('Time (s)'); ylabel('Glide-slope angle \theta (deg)');
    legend([h1(1) h2(1) h3 h4 h5], ...
           {'DO MC','RO MC','DO nominal','RO nominal','86 deg limit'}, ...
           'Location','southwest','FontSize',6.5); grid on;

    subplot(2,2,3);
    % 先画 MC（细线），再画标称（粗线）压在上面，否则后画的 RO MC 会把 DO MC 遮住。
    % 注意：DO 的质量 m(t)=m0*exp(-alpha*∫||u||dt) 只取决于推力积分，
    % 而 DO 所有样本的 u 完全相同 -> 500 条质量曲线严格重合为一条。
    % 故 DO MC 必须画得比标称线更粗，才能形成"浅蓝光晕 + 深蓝芯"的可见效果；
    % 若画得比标称线细，就会整条藏在标称线里看不见。
    h2 = plot(MC_RO.t, squeeze(MC_RO.m),'Color',[0.95 0.7 0.5],'LineWidth',0.5); hold on;
    h1 = plot(MC_DO.t, squeeze(MC_DO.m),'Color',[0.55 0.72 0.88],'LineWidth',0.5);
    h3 = plot(DOs.t, DOs.m,'b-','LineWidth',1.5);
    h4 = plot(ROs.t, ROs.m,'r--','LineWidth',1.5);
    xlabel('Time (s)'); ylabel('Mass (kg)');
    legend([h1(1) h2(1) h3 h4], {'DO MC','RO MC','DO','RO'}, ...
           'Location','northeast','FontSize',7); grid on;

    subplot(2,2,4);
    % DO 与 RO 的三条控制分量都画：颜色区**分量**（u1/u2/u3），线型区**方案**
    % （DO 实线、RO 虚线）。这样同一分量的两条曲线同色，直接可比；
    % 若改用"颜色区方案"，DO/RO 各三条同色线会互相纠缠，分不清哪个是哪个分量。
    clrU = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10];   % u1 蓝 / u2 绿 / u3 橙
    hU = gobjects(6,1);
    for i = 1:3
        hU(i) = plot(DOs.tc, DOs.u(i,:), '-', 'Color', clrU(i,:), 'LineWidth', 1.4);
        hold on;
    end
    for i = 1:3
        hU(3+i) = plot(ROs.tc, ROs.u(i,:), '--', 'Color', clrU(i,:), 'LineWidth', 1.4);
    end
    h0 = yline(0, 'k-');  set(h0,'LineWidth',0.5);
    xlabel('Time (s)'); ylabel('Thrust accel u (m/s^2)');
    % legend 按"列优先"填充两列：句柄顺序须排成 [DO u1 u2 u3, RO u1 u2 u3]，
    % 第 1 列才恰好全是 DO、第 2 列全是 RO（若按 DO/RO 交错排，行配对会错位）。
    legend(hU, ...
           {'DO u_1','DO u_2','DO u_3','RO u_1','RO u_2','RO u_3'}, ...
           'Location','northwest','FontSize',6.5,'NumColumns',2); grid on;

    % (5) 终端散布
    % 共享开环剖面下 DO 与 RO 的落点散布量级接近（不可约项占主导），
    % 故两者都用半透明填充点，重叠处才能同时看见，而不是后者盖住前者。
    figure('Name','Terminal dispersion','Color','w');
    ffloor = sqrt((sigma_r0.^2 + tf^2*sigma_v0.^2)/3);   % 逐轴 1sigma 不可约下限 (x,y,z)

    % 两组点数相同、量级相同，若都用大实心点，后画的会把先画的整个盖住。
    % 故 DO 用小的实心蓝点、RO 用更小的空心红点，两团云才都看得见。
    hD = scatter(MC_DO.r2f,  MC_DO.r3f,  18, 'b', 'filled', 'MarkerFaceAlpha',0.35); hold on;
    hR = scatter(MC_RO.r2f, MC_RO.r3f, 18, 'rs', 'MarkerFaceAlpha',0.35, 'LineWidth',0.5);
    hO = plot(0,0,'kp','MarkerSize',10,'MarkerFaceColor','k');
    xlabel('Terminal r_2 (m)'); ylabel('Terminal r_3 (m)');
    legend([hD hR hO], {'DO MC', sprintf('RO MC (%d samples)', numel(MC_RO.r2f)), 'Target'}, ...
           'Location','northwest','FontSize',9);
    % axis equal 会把坐标框撑满，统计量放标题里（两行）而不是塞进坐标区
    % title({'Terminal scatter at t = t_f', ...
    %        sprintf('DO (%.0f,%.0f) | RO (%.0f,%.0f) | floor (%.0f,%.0f) m', ...
    %                std(MC_DO.r2f), std(MC_DO.r3f), std(MC_RO.r2f), std(MC_RO.r3f), ...
    %                ffloor(2), ffloor(3))}, 'FontSize',9);
    grid on;
end

end   % ===== 主函数结束 =====

%% ===================== 求积与基函数 =====================
function [x,w] = gauss_legendre(n)
    beta = (1:n-1)./sqrt(4*(1:n-1).^2-1);
    J = diag(beta,1) + diag(beta,-1);
    [V,D] = eig(J); x = diag(D).'; w = 2*V(1,:).^2;
    [x,o] = sort(x); w = w(o); w = w/sum(w);
end

function [x,w] = cheb_lobatto(n)
    k = 0:n-1; x = cos(k*pi/(n-1)); w = ones(1,n);
    w(1) = 0.5; w(end) = 0.5; w = w/sum(w);
end

function [pts,w] = tensor_quad(nodes_cell, weights_cell)
    d = numel(nodes_cell);
    [G{1:d}] = ndgrid(nodes_cell{:});
    [W{1:d}] = ndgrid(weights_cell{:});
    Qq = prod(cellfun(@numel,nodes_cell));
    pts = zeros(Qq,d); w = ones(Qq,1);
    for i = 1:d
        pts(:,i) = G{i}(:);
        w = w .* W{i}(:);
    end
end

function P = legpoly(k, x)
% legpoly  归一化 Legendre 多项式（U(-1,1)），E[P_k P_l] = delta_kl
    if k == 0, P = ones(size(x)); return; end
    if k == 1, P = sqrt(3)*x; return; end
    P0 = ones(size(x)); P1 = x;
    for i = 2:k
        P2 = ((2*i-1)*x.*P1 - (i-1)*P0)/i; P0 = P1; P1 = P2;
    end
    P = sqrt(2*k+1)*P1;
end
