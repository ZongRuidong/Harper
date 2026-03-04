clc; clear; close all;

%% ========== 全局参数 ==========
d = 0.5; % 阵元间距(λ)
L = 2;   % 双源
sep_deg = 15; % 源间隔（近邻）
theta_true = [-sep_deg/2, +sep_deg/2];
theta_scan = -90:1:90; % 扫描栅格
SNR_dB = 10; % 固定 0 dB
SNR = 10^(SNR_dB/10);
M_fix = 16; % 固定阵元数
N_fix = 100; % 固定快拍数

%% ========== Fig2：RMSE vs SNR ==========
SNR_list = [-10 -5 0 10 20]; % SNR列表
curvesSNR = rmse_curves_vs_snr(SNR_list);

figure('Color','w','Name','Fig2 RMSE vs SNR'); hold on; box on; grid on;
plot(SNR_list, curvesSNR.MUSIC, '-o', 'LineWidth', 1.8, 'Color', [0 0 1]);       % MUSIC (Blue)
plot(SNR_list, curvesSNR.CAPON, '-s', 'LineWidth', 1.8, 'Color', [1 0.5 0]);     % Capon (Orange)
plot(SNR_list, curvesSNR.L1SVD, '-^', 'LineWidth', 1.8, 'Color', [0.5 0 0.5]);   % L1-SVD (Purple)
plot(SNR_list, curvesSNR.WLJ,   '-d', 'LineWidth', 1.8, 'Color', [1 1 0]);       % WLJ-ADMM (Yellow)
xlabel('SNR (dB)');
ylabel('RMSE (deg)');
title(sprintf('Fig2 RMSE vs SNR (M=%d, N=%d, sep=%d^\\circ)', M_fix, N_fix, sep_deg));
legend('MUSIC', 'Capon', 'L1-SVD', 'WLJ-ADMM', 'Location', 'northeast');
saveas(gcf, 'Fig2 RMSE vs SNR.png');
%% ========== Fig3：RMSE vs sep ==========
sep_list = [2 5 10 15 20]; % sep列表
curvesSEP = rmse_curves_vs_sep(sep_list); 

figure('Color','w','Name','Fig3 RMSE vs sep'); hold on; box on; grid on;
plot(sep_list, curvesSEP.MUSIC, '-o', 'LineWidth', 1.8, 'Color', [0 0 1]);       % MUSIC (Blue)
plot(sep_list, curvesSEP.CAPON, '-s', 'LineWidth', 1.8, 'Color', [1 0.5 0]);     % Capon (Orange)
plot(sep_list, curvesSEP.L1SVD, '-^', 'LineWidth', 1.8, 'Color', [0.5 0 0.5]);   % L1-SVD (Purple)
plot(sep_list, curvesSEP.WLJ,   '-d', 'LineWidth', 1.8, 'Color', [1 1 0]);       % WLJ-ADMM (Yellow)
xlabel('Separation (deg)');
ylabel('RMSE (deg)');
title(sprintf('Fig3 RMSE vs Separation (M=%d, N=%d, SNR=%d dB)', M_fix, N_fix, SNR_dB));
legend('MUSIC', 'Capon', 'L1-SVD', 'WLJ-ADMM', 'Location', 'northeast');
saveas(gcf, 'Fig3 RMSE vs Separation.png');
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

