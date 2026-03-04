clear; clc; close all;

%% ======== 参数设置 =========
M = 16; N = 500; d = 0.5;
theta_true = [-5, 5];
theta_scan = -90:0.5:90;
K_scan = length(theta_scan); % 防止变量名冲突
SNR_dB = -5;
SNR = 10^(SNR_dB / 10);

%% ======== 构造导向矩阵 (用户原始逻辑 0:M-1) =========
A = zeros(M, K_scan);
for k = 1:K_scan
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

%% ==== 调用 MUSIC (已有) ====
opts = struct('diagload', 1e-3, 'normalize', true, 'remove_mean', true, ...
              'minpkdist_deg', 2, 'prominence', [], 'smooth_win', 1);
% 注意：假设你已有 music_doa 函数，如果没有需自行补充，这里仅保留调用
try
    [P_music, ~] = music_doa(Y, M, 2, theta_scan, A, opts);
catch
    P_music = zeros(size(theta_scan)); % 占位
end

%% ==== 调用 Capon (已有) ====
R = (Y * Y') / N;
try
    P_capon = capon_doa(R, A);
    P_capon = P_capon / max(P_capon);
catch
    P_capon = zeros(size(theta_scan));
end

%% ==== 调用 L1-SVD (已有) ====
try
    params_l1 = struct('lambda1', 0.02, 'keep_energy', 0.95, ...
                       'rho', 1, 'iters', 30, 'L', 2);
    [~, P_l1] = l1svd_doa(Y, theta_scan, A, params_l1);
    P_l1 = P_l1 / max(P_l1);
catch
    P_l1 = zeros(size(theta_scan));
end

%% ==== 调用 WLJ-ADMM (已有) ====
try
    params_wlj = struct('lambda1', 0.01, 'lambdaJ', 0.002, 'epsJ', 1e-6, ...
                        'epsL', 1e-6, 'r_energy', 0.95, ...
                        'rho', 1, 'iters', 30, 'mm_outer', 2, 'L', 2);
    [~, P_wlj] = wljadmm_doa(Y, theta_scan, A, R, params_wlj);
    P_wlj = P_wlj / max(P_wlj);
catch
    P_wlj = zeros(size(theta_scan));
end

%% ==== 【新增】调用 RV-MCP-L1-SVD (方案三) ====
% 直接调用封装好的子函数
fprintf('Running RV-MCP-L1-SVD...\n');
[~, P_rv_mcp] = DOA_RV_MCP_L1_SVD(Y, 2, theta_scan, d);

%% ==== 绘图与可视化 ====
% 归一化小工具
norm01  = @(x) x ./ max(x(:) + eps);
adjust  = @(P,cfg) norm01( max( cfg.gain * (norm01(P).^cfg.gamma), cfg.floor ) );

% --- 可视化参数微调 ---
viz.l1.gain = 1.0; viz.l1.gamma = 0.80; viz.l1.floor = 0.0; viz.l1.smooth = 9;
viz.wlj.gain = 1.0; viz.wlj.gamma = 1.05; viz.wlj.floor = 0.0; viz.wlj.smooth = 5;

% 方案三的可视化参数 (MCP通常很稀疏，不需要太强的Gamma压扩)
viz.rv.gain   = 1.0; 
viz.rv.gamma  = 0.6;   % 稍微压一下让峰宽一点点，方便看
viz.rv.floor  = 0.0; 
viz.rv.smooth = 3;     % 轻微平滑

% 应用调整
P_l1  = adjust(P_l1 , viz.l1 );
P_wlj = adjust(P_wlj, viz.wlj);
P_rv_mcp = adjust(P_rv_mcp, viz.rv);

% 平滑处理
if viz.l1.smooth  > 1, P_l1  = movmean(P_l1 , viz.l1.smooth ); end
if viz.wlj.smooth > 1, P_wlj = movmean(P_wlj, viz.wlj.smooth); end
if viz.rv.smooth  > 1, P_rv_mcp = movmean(P_rv_mcp, viz.rv.smooth); end

% 最终归一化
P_l1  = P_l1  ./ max(P_l1  + eps);
P_wlj = P_wlj ./ max(P_wlj + eps);
P_rv_mcp = P_rv_mcp ./ max(P_rv_mcp + eps);

%% 绘图
figure('Color','w', 'Position', [100, 100, 800, 500]);
plot(theta_scan, P_music, 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1.0); hold on; % MUSIC 灰色垫底
plot(theta_scan, P_capon, 'o:', 'LineWidth', 1.2);
plot(theta_scan, P_l1, 'b-.', 'LineWidth', 1.5);
plot(theta_scan, P_wlj, 'k:', 'LineWidth', 1.8);
plot(theta_scan, P_rv_mcp, 'r-', 'LineWidth', 2); % 方案三 红色实线突出

% 真实 DOA 虚线
for th = theta_true
    xline(th, 'k--', 'LineWidth', 1.2);
end

legend('MUSIC','CAPON','L1-SVD','WLJ-ADMM','RV-MCP-L1-SVD','True DOA','Location','northeast');
xlabel('\theta (deg)'); ylabel('Normalized Spatial Spectrum');
title(['DOA Comparison (Coherent, SNR = ' num2str(SNR_dB) ' dB)']);
xlim([-20 20]); % 聚焦看中心区域
grid on;
