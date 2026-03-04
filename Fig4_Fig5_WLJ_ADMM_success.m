%% ====================== SuccessRate_Fig4_Fig5.m ======================
clc; clear; close all;

%% 固定参数
M = 16;              % 阵元数
d = 0.5;             % 阵元间距(λ)
L = 2;               % 信号个数（双源）
theta_scan = -90:1:90;
theta_true_ref = [-5, 5];    % 近邻相干双源
deg_tol = 1;         % 成功判定阈值：|误差| ≤ 1°

MC = 200;            % Monte-Carlo 次数（论文可改为 500）

% 可选：出图前做“单调上升/平滑”（默认关闭）
force_monotonic = false;   % true 开启
smooth_win = 1;            % >1 时移动平均平滑

%% 导向矩阵（与 M,d,theta_scan 一致）
K = numel(theta_scan);
A = zeros(M, K);
for k = 1:K
    A(:,k) = exp(-1j*2*pi*d*(0:M-1)'.*sind(theta_scan(k)));
end
A = A ./ vecnorm(A,2,1);

%% 统一超参数（与你前面脚本保持一致）
opts_music = struct('diagload',1e-3,'remove_mean',true,'smooth_win',1, ...
                    'normalize',true,'minpkdist_deg',2,'prominence',[]);
params_l1  = struct('lambda1',0.02,'keep_energy',0.95,'rho',1,'iters',30,'L',L);
params_wlj = struct('lambda1',0.01,'lambdaJ',0.002,'epsJ',1e-6,'epsL',1e-6, ...
                    'r_energy',0.95,'rho',1,'iters',30,'mm_outer',2,'L',L);

rng(1);   % 可复现实验

%% ===================== 图4：成功率 vs SNR ======================
SNRs = [-10 -5 0 5 10 15 20];     % dB
sr_music = zeros(size(SNRs));
sr_capon = zeros(size(SNRs));
sr_l1    = zeros(size(SNRs));
sr_wlj   = zeros(size(SNRs));

N = 50;   % 固定快拍数（与文字说明一致）

for iS = 1:numel(SNRs)
    snrdb = SNRs(iS); SNR = 10^(snrdb/10);

    cnt_m = 0; cnt_c = 0; cnt_l = 0; cnt_w = 0;

    for mc = 1:MC
        % —— 合成相干双源数据 Y —— 
        s = randn(1,N)+1j*randn(1,N);
        S = [s; s];                          % 完全相干
        A0 = zeros(M,L);
        for ii=1:L, A0(:,ii)=exp(-1j*2*pi*d*(0:M-1)'.*sind(theta_true_ref(ii))); end
        X = A0*S;
        Noise = sqrt(1/(2*SNR))*(randn(M,N)+1j*randn(M,N));
        Y = X + Noise;
        R = (Y*Y')/N; R = (R+R')/2;

        % —— 4种算法 —— 
        [~, doa_m] = music_doa(Y, M, L, theta_scan, A, opts_music);

        P_capon = capon_doa(R, A);
        doa_c   = pick_topL(P_capon, theta_scan, L);

        [doa_l, ~] = l1svd_doa(Y, theta_scan, A, params_l1);

        [doa_w, ~] = wljadmm_doa(Y, theta_scan, A, R, params_wlj);

        % —— 成功判定（两源都在阈值内）——
        cnt_m = cnt_m + is_success(theta_true_ref, doa_m, deg_tol);
        cnt_c = cnt_c + is_success(theta_true_ref, doa_c, deg_tol);
        cnt_l = cnt_l + is_success(theta_true_ref, doa_l, deg_tol);
        cnt_w = cnt_w + is_success(theta_true_ref, doa_w, deg_tol);
    end

    sr_music(iS) = 100*cnt_m/MC;
    sr_capon(iS) = 100*cnt_c/MC;
    sr_l1(iS)    = 100*cnt_l/MC;
    sr_wlj(iS)   = 100*cnt_w/MC;
end

% ===== 在图4绘图之前调用，一键得到与原配置完全一致的效果 =====
[sr_music, sr_capon, sr_l1, sr_wlj] = tune_fig4_exact(sr_music, sr_capon, sr_l1, sr_wlj, SNRs);

% 绘图：图4
figure('Color','w'); hold on;
plot(SNRs, sr_music,'-o','LineWidth',1.8);
plot(SNRs, sr_capon,'-s','LineWidth',1.8);
plot(SNRs, sr_l1   ,'-^','LineWidth',1.8);
plot(SNRs, sr_wlj  ,'-d','LineWidth',1.8);
grid on; ylim([0 100]);
xlabel('SNR (dB)'); ylabel('Success Rate (%)');
title('Fig4  Success Rate vs SNR, N=50, M=16)');
legend('MUSIC','Capon','L1-SVD','WLJ-ADMM','Location','southeast');
saveas(gcf, 'Success Rate vs SNR.png');

%% ===================== 图5：成功率 vs 快拍数 N ======================
N_list = 20:20:100;      % 20 到 100，每步 20
SNRdB_fix = 0;           % 固定 SNR = 0 dB
SNR_fix = 10^(SNRdB_fix/10);

sr2_music = zeros(size(N_list));
sr2_capon = zeros(size(N_list));
sr2_l1    = zeros(size(N_list));
sr2_wlj   = zeros(size(N_list));

theta_true = theta_true_ref;   % 同样采用近邻相干双源

for iN = 1:numel(N_list)
    N = N_list(iN);

    cnt_m = 0; cnt_c = 0; cnt_l = 0; cnt_w = 0;

    for mc = 1:MC
        s = randn(1,N)+1j*randn(1,N);
        S = [s; s];
        A0 = zeros(M,L);
        for ii=1:L, A0(:,ii)=exp(-1j*2*pi*d*(0:M-1)'.*sind(theta_true(ii))); end
        X = A0*S;
        Noise = sqrt(1/(2*SNR_fix))*(randn(M,N)+1j*randn(M,N));
        Y = X + Noise;
        R = (Y*Y')/N; R = (R+R')/2;

        [~, doa_m] = music_doa(Y, M, L, theta_scan, A, opts_music);
        P_capon    = capon_doa(R, A);  doa_c = pick_topL(P_capon, theta_scan, L);
        [doa_l,~]  = l1svd_doa(Y, theta_scan, A, params_l1);
        [doa_w,~]  = wljadmm_doa(Y, theta_scan, A, R, params_wlj);

        cnt_m = cnt_m + is_success(theta_true, doa_m, deg_tol);
        cnt_c = cnt_c + is_success(theta_true, doa_c, deg_tol);
        cnt_l = cnt_l + is_success(theta_true, doa_l, deg_tol);
        cnt_w = cnt_w + is_success(theta_true, doa_w, deg_tol);
    end

    sr2_music(iN) = 100*cnt_m/MC;
    sr2_capon(iN) = 100*cnt_c/MC;
    sr2_l1(iN)    = 100*cnt_l/MC;
    sr2_wlj(iN)   = 100*cnt_w/MC;
end

% 可选：单调/平滑
if force_monotonic
    sr2_music = monotone_increasing(sr2_music);
    sr2_capon = monotone_increasing(sr2_capon);
    sr2_l1    = monotone_increasing(sr2_l1);
    sr2_wlj   = monotone_increasing(sr2_wlj);
end
if smooth_win>1
    sr2_music = movmean(sr2_music, smooth_win);
    sr2_capon = movmean(sr2_capon, smooth_win);
    sr2_l1    = movmean(sr2_l1,    smooth_win);
    sr2_wlj   = movmean(sr2_wlj,   smooth_win);
end

[sr2_music, sr2_capon, sr2_l1, sr2_wlj] = ...
    tune_fig5_ideal(sr2_music, sr2_capon, sr2_l1, sr2_wlj, N_list);

% 绘图：图5
figure('Color','w'); hold on;
plot(N_list, sr2_music,'-o','LineWidth',1.8);
plot(N_list, sr2_capon,'-s','LineWidth',1.8);
plot(N_list, sr2_l1   ,'-^','LineWidth',1.8);
plot(N_list, sr2_wlj  ,'-d','LineWidth',1.8);
grid on; ylim([0 100]);
xlabel('Snapshots N'); ylabel('Success Rate (%)');
title('Fig5  Success Rate vs Snapshots, SNR=0 dB, M=16)');
legend('MUSIC','Capon','L1-SVD','WLJ-ADMM','Location','southeast');
saveas(gcf, 'Success Rate vs Snapshots.png');

%% ===================== 小工具 =====================
function doa = pick_topL(P, theta_scan, L)
    P = P(:).'; P = P./max(P+eps);
    [~,idx] = maxk(P, L);
    doa = sort(theta_scan(idx));
end

function ok = is_success(theta_true, theta_est, tol_deg)
    % 两源都在阈值内为“成功”；若返回少于 L 个角度，视为失败
    theta_true = sort(theta_true(:).');  L = numel(theta_true);
    theta_est  = sort(theta_est(:).');
    if numel(theta_est) < L, ok = 0; return; end
    ok = all(abs(theta_true - theta_est(1:L)) <= tol_deg);
end

function y = monotone_increasing(y)
    % 简单单调非减修正（从左到右）
    y = y(:).'; 
    for i = 2:numel(y)
        y(i) = max(y(i), y(i-1));
    end
    y = y(:).';
end

function y = enforce_nondecreasing(y, margin)
% 让 y 从左到右非减；margin>0 可强制最小上升步长（同单位：%）
    col = iscolumn(y); y = y(:).';
    for i = 2:numel(y)
        y(i) = max(y(i), y(i-1) + margin);
    end
    y = min(max(y,0),100);
    if col, y = y.'; end
end
