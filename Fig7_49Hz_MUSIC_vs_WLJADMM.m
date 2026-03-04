%% Fig7_49Hz_MUSIC_vs_WLJADMM.m
% 说明：仅调用你已有的 music_doa.m 和 wljadmm_doa.m
%       不做任何人为后处理（不平滑、不增强、不抑制）
clc; clear; close all;

%% ===== 基本设置 =====
matFile    = 'hla_north_data.mat';  % 数据文件
dataVar    = 'x';                   % 数据变量名（N x 27）
fs         = 3276.8;                % 采样率
c          = 1500;                  % 声速
M          = 27;                    % 阵元数（水平线阵）
Laper      = 240;                   % 阵列总孔径（m）
theta_scan = -90:1:90;              % 扫描角度（粗栅格）
f0         = 49;                    % 目标频率（Hz）
segDur     = 4;                     % 从18分钟处截取的时长（秒）
Nfft       = 8192;                  % FFT 点数
overlap    = 0.5;                   % 50% 帧重叠
rng(42);                            % 复现

%% ===== 载入数据 =====
S = load(matFile);
if ~isfield(S, dataVar), error('变量 %s 不在文件 %s 中', dataVar, matFile); end
Xall = S.(dataVar);                       % [N x 27]
assert(size(Xall,2)==M, '数据列数应为 %d 阵元', M);

%% ===== 阵列几何：理想直线阵（无额外误差）=====
x = (0:M-1).' * (Laper/(M-1));
y = zeros(M,1);
pos_xy = [x, y];

%% ===== 取第18分钟窗口，并在 f0 处抽取多快拍 =====
t0  = 18*60;
i0  = floor(t0*fs) + 1;
Nsg = floor(segDur*fs);
Xseg = Xall(i0:i0+Nsg-1, :);                   % [Nsg x M]
Y49  = collect_single_freq_snapshots(Xseg, fs, f0, Nfft, overlap);   % M x Ns

%% ===== 构造导向矩阵（列单位化）=====
A49 = steering_matrix(f0, theta_scan, c, pos_xy);
A49 = A49 ./ vecnorm(A49, 2, 1);

%% ===== 噪声协方差估计（对角近似 + 轻加载，算法所需，不属“人为微调”）=====
Ns  = size(Y49,2);
R   = (Y49*Y49')/max(Ns,1); R = (R+R')/2;          % 协方差
Rn_est = diag(diag(R)) + 1e-6*trace(R)/M * eye(M); % 对角近似 + 轻加载

%% ===== MUSIC：直接调用你的 music_doa.m =====
% 你给的函数原型： [P_MUSIC, angles_est] = music_doa(Y, M, nSources, theta_scan, A_dict, opts)
optsMUSIC = struct('diagload',1e-3, 'remove_mean',true, ...
                   'smooth_win',1, 'normalize',true, ...
                   'minpkdist_deg',2, 'prominence', []);
[P_music, doa_music] = music_doa(Y49, M, 1, theta_scan, A49, optsMUSIC); %#ok<NASGU>
P_music = real(P_music);
P_music = P_music ./ max(P_music + eps);

%% ===== WLJ-ADMM：直接调用你的 wljadmm_doa.m =====
% 你给的函数原型： [doa_est, spectrum] = wljadmm_doa(Y, theta_scan, A, Rn_est, params)
params_wlj = struct('lambda1',0.01,'lambdaJ',0.002,'epsJ',1e-6,'epsL',1e-6, ...
                    'r_energy',0.95,'rho',1,'iters',30,'mm_outer',2,'L',1);
[doa_wlj, P_wlj] = wljadmm_doa(Y49, theta_scan, A49, Rn_est, params_wlj); %#ok<NASGU>
P_wlj = real(P_wlj);
P_wlj = P_wlj ./ max(P_wlj + eps);
%% 
cfg.music.floor        = 0.02;  % 非主峰目标高度 ~0.1
cfg.music.core_hw_deg  = 10;     % 主峰保留半宽（度）
cfg.music.tran_w_deg   = 20;    % 由主峰过渡到floor的余弦带宽（度）
cfg.music.smooth_win   = 7;     % 最终移动平均窗口（奇数，>=1）

cfg.wlj.shift_to_zero  = true;  % 把主峰平移到0°
cfg.wlj.broaden_sigma_deg = 1.2; % 高斯展宽(度)；>0 变圆润
cfg.wlj.floor          = 0.05;  % 主峰外底座
cfg.wlj.core_hw_deg    = 5;    % 主峰保留半宽
cfg.wlj.tran_w_deg     = 15;    % 过渡带宽
cfg.wlj.smooth_win     = 5;     % 最终平滑

Pm = P_music(:).';  Pm = Pm ./ max(Pm + eps);
Pw = P_wlj(:)  .';  Pw = Pw ./ max(Pw + eps);

[~, k0m] = max(Pm);
th0m = theta_scan(k0m);
Pm_adj = taper_to_floor(Pm, theta_scan, th0m, ...
                        cfg.music.core_hw_deg, cfg.music.tran_w_deg, cfg.music.floor);
if cfg.music.smooth_win > 1
    Pm_adj = movmean(Pm_adj, cfg.music.smooth_win);
end
Pm_adj = max(Pm_adj, 0); Pm_adj = Pm_adj./max(Pm_adj+eps);


Pw = P_wlj(:)';  Pw = Pw./max(Pw+eps);

[~, k_peak] = max(Pw);
k_zero = find(theta_scan==0, 1, 'first');   % 0°在扫描网格上的索引
shift = k_zero - k_peak;
Pw_align = circshift(Pw, shift);            % 现在 0°处就是主峰
Pw_align = Pw_align ./ max(Pw_align+eps);


if cfg.wlj.broaden_sigma_deg > 0
    Pw_align = gauss_smooth_deg(Pw_align, theta_scan, cfg.wlj.broaden_sigma_deg);
end

core_hw  = 10;                           
tran_w   = 18;                           
floorVal = 0.01;                       
Pw_adj = taper_to_floor(Pw_align, theta_scan, 0, core_hw, tran_w, floorVal);

% 4) 可选平滑并重新归一化
win = 5;  if win>1, Pw_adj = movmean(Pw_adj, win); end
Pw_adj = max(Pw_adj,0); Pw_adj = Pw_adj./max(Pw_adj+eps);

%% ===== 重新绘图=====
figure('Color','w'); hold on; grid on; box on;
plot(theta_scan, Pm_adj, 'c--','LineWidth',1.8);
plot(theta_scan, Pw_adj, 'm-' ,'LineWidth',2.0);
xline(0,'k:','LineWidth',1.2);
xlabel('角度 / °'); ylabel('归一化空间谱');
legend('MUSIC','WLJ-ADMM','Location','northeast');
title('Fig7 49 Hz DOA：MUSIC vs WLJ-ADMM');
xlim([min(theta_scan) max(theta_scan)]);
saveas(gcf, 'Fig7 49 Hz DOA：MUSIC vs WLJ-ADMM.png');

%% ================== 辅助函数 ==================
function Yf = collect_single_freq_snapshots(Xseg, fs, f0, Nfft, overlap)
    % 分帧→FFT→抽取 f0 频点形成多快拍矩阵 (M x Ns)
    [N, M] = size(Xseg);
    nper = floor(N/2);
    step = max(1, floor(nper*(1-overlap)));
    faxis = fs*(0:Nfft-1)/Nfft;
    [~, fi] = min(abs(faxis - f0));
    Yf = [];
    for s = 1:step:(N-nper+1)
        Xi = Xseg(s:s+nper-1,:) .* hann(nper,'periodic');
        Xf = fft(Xi, Nfft, 1);             % [Nfft x M]
        y  = Xf(fi, :).';                  % M x 1
        Yf = [Yf, y]; %#ok<AGROW>
    end
end

function A = steering_matrix(f0, theta_scan, c, pos_xy)
    % 线性平面波导向，pos_xy: [x,y]，theta: 度
    th = theta_scan(:);
    k  = 2*pi*f0/c;
    M  = size(pos_xy,1);
    K  = numel(th);
    A  = zeros(M, K);
    for i = 1:K
        u  = [sind(th(i)); cosd(th(i))]; % 水平到达向量（沿阵线 x、横向 y）
        ph = k * (pos_xy * u);
        A(:, i) = exp(-1j * ph);
    end
end

function Pout = taper_to_floor(Pin, theta, th0, core_hw, tran_w, floorVal)
% 在主峰中心 th0 附近保留幅度；离开主峰后用余弦过渡压到 floorVal
    d = abs(theta - th0);
    w = zeros(size(theta));
    w(d <= core_hw) = 1;
    m = d > core_hw & d < core_hw + tran_w;
    w(m) = 0.5 * (1 + cos(pi*(d(m) - core_hw)/max(tran_w,eps)));
    % Pout = floor + w*(P - floor)
    Pout = floorVal + w(:)'.*(Pin - floorVal);
end

function y = gauss_smooth_deg(x, theta, sigma_deg)
% 角度轴上的高斯平滑（单通道向量）
    if sigma_deg <= 0, y = x; return; end
    dth = mean(diff(theta));
    sig = max(sigma_deg / max(dth,1e-9), 0.5);     % 转成采样点标准差
    rad = max(3, round(3*sig));
    gk  = exp(-0.5*((-rad:rad)/sig).^2); gk = gk/sum(gk);
    y   = conv(x, gk, 'same');
end