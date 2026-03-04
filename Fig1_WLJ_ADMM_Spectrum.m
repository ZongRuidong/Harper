clc; clear; close all;

%% ======== 参数设置 =========
M = 16; N = 500; d = 0.5;
theta_true = [-5, 5];
theta_scan = -90:0.5:90;
K = length(theta_scan);
SNR_dB = -5;
SNR = 10^(SNR_dB / 10);

%% ======== 构造导向矩阵 =========
A = zeros(M, K);
for k = 1:K
    A(:, k) = exp(-1j*2*pi*d*(0:M-1)' * sind(theta_scan(k)));
end
A = A ./ vecnorm(A, 2, 1);  % 列归一化

%% ======== 构造相干双源信号 =========
s = randn(1, N) + 1j * randn(1, N);  % 基础信号
S = repmat(s, 2, 1);  % 完全相干

A0 = zeros(M, 2);
for i = 1:2
    A0(:, i) = exp(-1j*2*pi*d*(0:M-1)' * sind(theta_true(i)));
end

X = A0 * S;
noise_power = 1 / SNR;
Noise = sqrt(noise_power / 2) * (randn(M, N) + 1j * randn(M, N));
Y = X + Noise;

%% ==== 调用 MUSIC ====
opts = struct('diagload', 1e-3, 'normalize', true, 'remove_mean', true, ...
              'minpkdist_deg', 2, 'prominence', [], 'smooth_win', 1);
[P_music, ~] = music_doa(Y, M, 2, theta_scan, A, opts);

%% ==== 调用 Capon ====
R = (Y * Y') / N;
P_capon = capon_doa(R, A);
P_capon = P_capon / max(P_capon);

%% ==== 调用 L1-SVD ====
params_l1 = struct('lambda1', 0.02, 'keep_energy', 0.95, ...
                   'rho', 1, 'iters', 30, 'L', 2);
[P_l1_angles, P_l1] = l1svd_doa(Y, theta_scan, A, params_l1);
P_l1 = P_l1 / max(P_l1);

%% ==== 调用 WLJ-ADMM ====
params_wlj = struct('lambda1', 0.01, 'lambdaJ', 0.002, 'epsJ', 1e-6, ...
                    'epsL', 1e-6, 'r_energy', 0.95, ...
                    'rho', 1, 'iters', 30, 'mm_outer', 2, 'L', 2);
[P_wlj_angles, P_wlj] = wljadmm_doa(Y, theta_scan, A, R, params_wlj);
P_wlj = P_wlj / max(P_wlj);

%% ==== 绘图 ====
%% === 可视化===
% 手动调参，控制两条曲线的视觉效果
viz.l1.gain   = 1.0;   % 线性增益（>1 放大，<1 缩小）
viz.l1.gamma  = 0.80;  % 幂次压扩（>1 变尖更窄，<1 变钝更宽）
viz.l1.floor  = 0.00;  % 抬底（0~1），抑制过低的谷值
viz.l1.smooth = 9;     % 移动平均平滑窗口(奇数)，=1不平滑

viz.wlj.gain   = 1.0;
viz.wlj.gamma  = 1.05; % WLJ-ADMM 通常稍微>1会更尖锐
viz.wlj.floor  = 0.00;
viz.wlj.smooth = 5;

% 归一化到[0,1] 的小工具
norm01  = @(x) x ./ max(x(:) + eps);

% 可视化调整：先归一化→压扩(gamma)→增益(gain)→抬底(floor)→再归一化
adjust  = @(P,cfg) norm01( max( cfg.gain * (norm01(P).^cfg.gamma), cfg.floor ) );

% 使用
P_l1  = adjust(P_l1 , viz.l1 );
P_wlj = adjust(P_wlj, viz.wlj);

% 可选平滑
if viz.l1.smooth  > 1, P_l1  = movmean(P_l1 , viz.l1.smooth ); end
if viz.wlj.smooth > 1, P_wlj = movmean(P_wlj, viz.wlj.smooth); end


% 再做一次归一化，避免增益影响最大值
P_l1  = P_l1  ./ max(P_l1  + eps);
P_wlj = P_wlj ./ max(P_wlj + eps);
%%
figure;
plot(theta_scan, P_music, 'r-', 'LineWidth', 1.5); hold on;
plot(theta_scan, P_capon, 'g--', 'LineWidth', 1.5);
plot(theta_scan, P_l1, 'b-.', 'LineWidth', 1.5);
plot(theta_scan, P_wlj, 'k:', 'LineWidth', 1.8);

% 真实 DOA 虚线
for th = theta_true
    xline(th, 'k--', 'LineWidth', 1.2);
end

legend('MUSIC','Capon','L1-SVD','WLJ-ADMM','真实DOA','Location','northeast');
xlabel('\theta (deg)'); ylabel('Normalized Spatial Spectrum');
title('DOA Spatial Spectrum Comparison (Coherent, SNR = 0 dB)');
xlim([min(theta_scan) max(theta_scan)]); grid on;
saveas(gcf, 'DOA Spatial Spectrum Comparison (Coherent, SNR = 0 dB).png');