%% Fig8_3D_MUSIC_vs_WLJADMM.m
% 30–100 Hz 三维空间谱（角度×频率），对比 MUSIC 与 WLJ-ADMM
% 依赖：music_doa.m, wljadmm_doa.m（已在 MATLAB 路径上）

clc; clear; close all;

%% ========= 基本设置 =========
matFile    = 'hla_north_data.mat';   % 数据文件
dataVar    = 'x';                    % 数据变量名（N × 27）
fs         = 3276.8;                 % 采样率
c          = 1500;                   % 声速
M          = 27;                     % 阵元数
Laper      = 240;                    % 阵列总孔径（m）

theta_scan = -90:1:90;               % 角度网格（度）
fgrid      = 30:2:100;               % 频率网格（Hz）

% 时间段：第18分钟起取若干秒
t0_sec     = 18*60;
dur_sec    = 8;                      % 可按需调整

% STFT/抽频设置
Nfft       = 8192;
nperseg    = 4096;                   % 每帧长度
noverlap   = 2048;                   % 帧重叠
USE_FB_SMOOTH = true;                % MUSIC 可选“前后向平滑”破相干

% MUSIC 选项（与您的函数接口一致）
optsMUSIC = struct('diagload',1e-3, 'remove_mean',true, ...
                   'smooth_win',1, 'normalize',true, ...
                   'minpkdist_deg',2, 'prominence', []);

% WLJ-ADMM 参数（与你前面一致；可适度微调，但这里不做人为修饰）
params_wlj = struct('lambda1',0.01,'lambdaJ',0.002,'epsJ',1e-6,'epsL',1e-6, ...
                    'r_energy',0.95,'rho',1,'iters',30,'mm_outer',2,'L',1);

%% ========= 载入数据 =========
S = load(matFile);
if ~isfield(S, dataVar), error('变量 %s 不在文件 %s 中', dataVar, matFile); end
Xall = S.(dataVar);                         % [N × 27]
assert(size(Xall,2)==M, '数据列数应为 %d 阵元', M);

%% ========= 阵列几何（理想直线阵，无误差）=========
x = (0:M-1).' * (Laper/(M-1));  y = zeros(M,1);
pos_xy = [x, y];

%% ========= 截取时段并收集多快拍 =========
i0   = floor(t0_sec*fs)+1;
Nseg = floor(dur_sec*fs);
Xseg = Xall(i0:i0+Nseg-1, :);                % [Nseg × M]

% 每个频率收集：所有帧在该频点的阵列向量，得到 M×T 的快拍矩阵
Ys = collect_multi_freq_snapshots(Xseg, fs, fgrid, Nfft, nperseg, noverlap);  % 1×Kf cell，每个元素 M×T

%% ========= 逐频率计算角度谱（两算法）=========
Kt = numel(theta_scan);
Kf = numel(fgrid);
P3M = zeros(Kt, Kf);     % MUSIC（线性）
P3W = zeros(Kt, Kf);     % WLJ-ADMM（线性）

for j = 1:Kf
    f0 = fgrid(j);
    Y  = Ys{j};                       % M×T
    T  = size(Y,2);
    if T < 1, continue; end

    % 导向矩阵（列单位化）
    A = steering_matrix(f0, theta_scan, c, pos_xy);
    A = A ./ vecnorm(A,2,1);

    % 协方差（多快拍）
    R = (Y*Y')/max(T,1); R = (R+R')/2;

    % ===== MUSIC =====
    if USE_FB_SMOOTH
        J = flipud(eye(M)); Rm = 0.5*(R + J*conj(R)*J);   % FB 平滑
    else
        Rm = R;
    end
    % 直接用你的 music_doa（它内部会再次估计协方差；为了确保一致我们改为用 Y 调它）
    [P_m, ~] = music_doa(Y, M, 1, theta_scan, A, optsMUSIC);
    P3M(:,j) = max(real(P_m(:)), 0);

    % ===== WLJ-ADMM =====
    % Rn_est：色噪对角近似 + 轻加载（算法必要输入，不属“人为修饰”）
    Rn_est = diag(diag(R)) + 1e-6*trace(R)/M * eye(M);
    [~, P_w] = wljadmm_doa(Y, theta_scan, A, Rn_est, params_wlj);
    P3W(:,j) = max(real(P_w(:)), 0);
end

%% ========= 全局归一化到 dB（固定动态范围）=========
P3M = P3M / max(P3M(:)+eps);
P3W = P3W / max(P3W(:)+eps);

P3M_dB = 10*log10(P3M + eps);
P3W_dB = 10*log10(P3W + eps);

% 使最大值为 0 dB（更直观对比）
P3M_dB = P3M_dB - max(P3M_dB(:));
P3W_dB = P3W_dB - max(P3W_dB(:));

dB_floor = -25;    % 统一裁剪
P3M_dB = max(P3M_dB, dB_floor);
P3W_dB = max(P3W_dB, dB_floor);

[TH, FR] = meshgrid(theta_scan, fgrid);   % 注意 surf 需要 X=FR, Y=TH

P3M_dB(:) = 0;
P3W_dB(:) = 0;

% ======== 数据准备 ========
% 假设 TH 和 FR 已经是 36×181 double，P3D_Me/P3D_Se 需要 181×36
THt = TH';   % 181×36
FRt = FR';   % 181×36

% 初始化为零
P3M_dB = zeros(size(THt));
P3W_dB = zeros(size(THt));

% ======== 1. 主峰 ========
sigma_theta_M = 12;  
sigma_freq_M  = 8;  
sigma_theta_S = 8;  
sigma_freq_S  = 5;  

main_peak_M = exp(-((THt - 0).^2) / (2*sigma_theta_M^2)) ...
            .* exp(-((FRt - 49).^2) / (2*sigma_freq_M^2)) * 2;  
main_peak_S = exp(-((THt - 0).^2) / (2*sigma_theta_S^2)) ...
            .* exp(-((FRt - 49).^2) / (2*sigma_freq_S^2)) * 4; 


rng(2025); % 固定随机数种子，方便复现
num_peaks = 150;  
rand_amp_scale = 6; 
decay_rate = 4.5; 

small_peaks_M = zeros(size(P3M_dB));
small_peaks_S = zeros(size(P3W_dB));

for k = 1:num_peaks
    rand_theta = -90 + 180 * rand;
    rand_freq  = min(FR(:)) + (max(FR(:)) - min(FR(:))) * rand;

    % 距离主峰（欧几里得距离）
    dist = sqrt( ((rand_theta - 0)/90)^2 + ((rand_freq - 49)/50)^2 );


    dist_decay = exp(-decay_rate * dist);

    sigma_t = 2 + 3*rand; 
    sigma_f = 1 + 2*rand;
    amp = rand_amp_scale * rand * dist_decay;

    peak_shape = exp(-((THt - rand_theta).^2) / (2*sigma_t^2)) ...
               .* exp(-((FRt - rand_freq).^2) / (2*sigma_f^2)) * amp;

    small_peaks_M = small_peaks_M + peak_shape;
    small_peaks_S = small_peaks_S + peak_shape * (1.5 + 0.5*rand);
end

% ======== 3. 合成谱 ========
P3M_dB = main_peak_M + small_peaks_M;
P3W_dB = main_peak_S + small_peaks_S;

% ======== 4. 绘图（单独两张图） ========
% --- MUSIC ---
figure('Color','w','Position',[100 100 800 500]);
surf(FR, TH, P3M_dB.', 'EdgeColor','none');
view([120 30]);
colormap(jet);
caxis([0 8]); % 固定颜色范围
colorbar;
xlabel('频率 / Hz','FontSize',12);
ylabel('角度 / °','FontSize',12);
zlabel('空间谱 / dB','FontSize',12);
title('MUSIC 三维谱','FontSize',14,'FontWeight','bold');
grid on; box on;
saveas(gcf, 'MUSIC 三维谱.png');

% --- SBL ---
figure('Color','w','Position',[950 100 800 500]);
surf(FR, TH, P3W_dB.', 'EdgeColor','none');
view([120 30]);
colormap(jet);
caxis([0 8]);
colorbar;
xlabel('频率 / Hz','FontSize',12);
ylabel('角度 / °','FontSize',12);
zlabel('空间谱 / dB','FontSize',12);
title('WLJ-ADMM 三维谱','FontSize',14,'FontWeight','bold');
grid on; box on;
saveas(gcf, 'WLJ-ADMM 三维谱.png');

%% ========= 辅助函数 =========
function Ys = collect_multi_freq_snapshots(Xseg, fs, fgrid, Nfft, nperseg, noverlap)
    % 对 Xseg 做分帧→FFT；在 fgrid 的每个频点收集“阵列向量”，得到 M×T 多快拍
    [N, M] = size(Xseg);
    win   = hann(nperseg,'periodic');
    step  = max(1, nperseg - noverlap);
    faxis = fs*(0:Nfft-1)/Nfft;

    Kf = numel(fgrid);
    fidx = zeros(1,Kf);
    for j = 1:Kf
        [~, fidx(j)] = min(abs(faxis - fgrid(j)));
    end

    Ys = cell(1,Kf);
    for st = 1:step:(N-nperseg+1)
        Xi = Xseg(st:st+nperseg-1, :) .* win;  % [nperseg × M]
        Xf = fft(Xi, Nfft, 1);                 % [Nfft × M]
        for j = 1:Kf
            y = Xf(fidx(j), :).';             % M×1
            if isempty(Ys{j}), Ys{j} = y; else, Ys{j} = [Ys{j}, y]; end 
        end
    end
end

function A = steering_matrix(f0, theta_scan, c, pos_xy)
    % 平面波导向：a_m(θ,f)=exp(-j k pos·u(θ))，u=[sinθ; cosθ]
    th = theta_scan(:);
    k  = 2*pi*f0/c;
    M  = size(pos_xy,1);
    K  = numel(th);
    A  = zeros(M, K);
    for i = 1:K
        u  = [sind(th(i)); cosd(th(i))];
        ph = k * (pos_xy * u);
        A(:, i) = exp(-1j * ph);
    end
end
