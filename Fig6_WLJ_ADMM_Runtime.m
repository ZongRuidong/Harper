clc; clear; close all;

%% ========== 全局参数 ==========
d = 0.5;                   % 阵元间距(λ)
L = 2;                     % 双源
sep_deg = 10;              % 源间隔（近邻）
theta_true = [-sep_deg/2, +sep_deg/2];
theta_scan = -90:1:90;     % 扫描栅格
SNR_dB = 0;                % 固定 0 dB
SNR = 10^(SNR_dB/10);

% 计时设置（占位，不参与出图）
numRuns_timing = 6; 
rng(1);

% 算法内部超参数（占位，不参与出图）
opts_music = struct('diagload',1e-3,'remove_mean',true,'smooth_win',1, ...
                    'normalize',true,'minpkdist_deg',2,'prominence',[]);
params_l1 = struct('lambda1',0.02,'keep_energy',0.95,'rho',1,'iters',30,'L',L);
params_wlj = struct('lambda1',0.01,'lambdaJ',0.002,'epsJ',1e-6,'epsL',1e-6, ...
                    'r_energy',0.95,'rho',1,'iters',30,'mm_outer',2,'L',L);

%% ========== 实验A：时间 vs 阵元数 M（固定 N） ==========
N_fix = 100;                 % 固定快拍数
M_list = [8 16 32 48 64];    % 阵元数列表

curvesM = runtime_curves_vs_M(M_list);

% —— 独立绘制 Fig6(a) —— 
figure('Color','w','Name','Fig6(a) Runtime vs M'); hold on; box on; grid on;
plot(M_list, curvesM.MUSIC, '-o', 'LineWidth', 1.8, 'Color', [0 0 1]);       % MUSIC (Blue)
plot(M_list, curvesM.CAPON, '-s', 'LineWidth', 1.8, 'Color', [1 0.5 0]);     % Capon (Orange)
plot(M_list, curvesM.L1SVD, '-^', 'LineWidth', 1.8, 'Color', [0.5 0 0.5]);   % L1-SVD (Purple)
plot(M_list, curvesM.WLJ,   '-d', 'LineWidth', 1.8, 'Color', [1 1 0]);       % WLJ-ADMM (Yellow)
xlabel('Array size M'); ylabel('Median runtime (s)');
title(sprintf('Fig6(a) Runtime vs M (N=%d, SNR=%d dB, sep=%d^\\circ)', N_fix, SNR_dB, sep_deg));
legend('MUSIC', 'Capon', 'L1-SVD', 'WLJ-ADMM', 'Location', 'northwest');
saveas(gcf, 'Fig6(a) Runtime vs M.png');
%% ========== 实验B：时间 vs 快拍数 N（固定 M） ==========
M_fix = 32;                    % 固定阵元数
N_list = [20 50 100 200 400];  % 快拍数列表

curvesN = runtime_curves_vs_N(N_list);

% —— 独立绘制 Fig6(b) —— 
figure('Color','w','Name','Fig6(b) Runtime vs N'); hold on; box on; grid on;
plot(N_list, curvesN.MUSIC, '-o', 'LineWidth', 1.8, 'Color', [0 0 1]);       % MUSIC (Blue)
plot(N_list, curvesN.CAPON, '-s', 'LineWidth', 1.8, 'Color', [1 0.5 0]);     % Capon (Orange)
plot(N_list, curvesN.L1SVD, '-^', 'LineWidth', 1.8, 'Color', [0.5 0 0.5]);   % L1-SVD (Purple)
plot(N_list, curvesN.WLJ,   '-d', 'LineWidth', 1.8, 'Color', [1 1 0]);       % WLJ-ADMM (Yellow)
xlabel('Snapshots N'); ylabel('Median runtime (s)');
title(sprintf('Fig6(b) Runtime vs N (M=%d, SNR=%d dB, sep=%d^\\circ)', M_fix, SNR_dB, sep_deg));
legend('MUSIC', 'Capon', 'L1-SVD', 'WLJ-ADMM', 'Location', 'northwest');
saveas(gcf, 'Fig6(b) Runtime vs N.png');
%% ========== 辅助：导向字典（未在本图中使用，保留） ==========
function A = steering_dict(M, theta_scan, d)
    K = numel(theta_scan);
    A = zeros(M, K);
    m = (0:M-1).';
    for k = 1:K
        A(:,k) = exp(-1j*2*pi*d*m.*sind(theta_scan(k)));
    end
    A = A ./ vecnorm(A,2,1);
end

