function [P_MUSIC, angles_est] = music_doa(Y, M, nSources, theta_scan, A_dict, opts)
% MUSIC_DOA  多信号分类(MUSIC)算法的DOA估计
%
% 用法：
%   [P_MUSIC, angles_est] = music_doa(Y, M, nSources, theta_scan, A_dict)
%   [P_MUSIC, angles_est] = music_doa(Y, M, nSources, theta_scan, A_dict, opts)
%
% 输入：
%   Y           - 观测数据 (M x N)，每列为一个快拍或频域快拍
%   M           - 阵元数（与 size(Y,1) 一致）
%   nSources    - 期望信源个数（用于构造噪声子空间与选峰）
%   theta_scan  - 扫描角度向量（度），如 -90:1:90
%   A_dict      - 导向矩阵 (M x K)，K = length(theta_scan)
%   opts        - 可选参数结构体（字段均可省略）
%       .diagload        对角加载系数（相对于 trace(R)/M 的比例），默认 1e-3
%       .remove_mean     是否去除通道均值（默认 true）
%       .smooth_win      谱平滑窗口（奇数，默认 1=不平滑）
%       .normalize       是否归一化到[0,1]（默认 true）
%       .minpkdist_deg   峰间最小角度间隔（度，默认 2）
%       .prominence      峰显著性（默认 [] 不约束）
%
% 输出：
%   P_MUSIC     - MUSIC 伪谱（1 x K）
%   angles_est  - 估计的 DOA 角度（1 x nSources）
%
% 说明：
%   该实现对协方差进行Hermitian化与对角加载；使用 1/(||E_n^H a||^2) 形式计算谱；
%   内置平滑与归一化，便于直接绘图；峰值用 findpeaks 自动选取 nSources 个。

    if nargin < 6, opts = struct(); end
    if ~isfield(opts, 'diagload'),      opts.diagload = 1e-3; end
    if ~isfield(opts, 'remove_mean'),    opts.remove_mean = true; end
    if ~isfield(opts, 'smooth_win'),     opts.smooth_win = 1; end
    if ~isfield(opts, 'normalize'),      opts.normalize = true; end
    if ~isfield(opts, 'minpkdist_deg'),  opts.minpkdist_deg = 2; end
    if ~isfield(opts, 'prominence'),     opts.prominence = []; end

    % 基本尺寸检查
    [Mm, N] = size(Y);
    if Mm ~= M
        error('music_doa: 输入 M 与 Y 的行数不一致。');
    end
    if size(A_dict,1) ~= M || size(A_dict,2) ~= numel(theta_scan)
        error('music_doa: A_dict 尺寸应为 M x length(theta_scan)。');
    end
    if nSources >= M
        error('music_doa: nSources 必须小于阵元数 M。');
    end

    % 1) 协方差矩阵估计
    Yw = Y;
    if opts.remove_mean
        Yw = Yw - mean(Yw, 2);        % 去除每通道直流分量
    end
    R = (Yw * Yw.') / max(N,1);       % M x M
    R = (R + R')/2;                   % Hermitian化

    % 2) 对角加载（提高稳健性）
    if opts.diagload > 0
        dl = opts.diagload * trace(R)/M;
        R = R + dl * eye(M);
    end

    % 3) 特征分解，构造噪声子空间
    [V, D] = eig(R);
    [~, ord] = sort(real(diag(D))), ord = ord(:); %#ok<NASGU>
    % 噪声子空间 = 最小 M-nSources 个特征值对应向量
    [~, idxAsc] = sort(real(diag(D)), 'ascend');
    En = V(:, idxAsc(1:M-nSources));      % M x (M-nSources)

    % 4) 计算 MUSIC 谱：P = 1 / ||E_n^H a||^2
    % 为避免显式循环，向量化实现
    EnH_A = En' * A_dict;                 % (M-nS) x K
    denom = sum(abs(EnH_A).^2, 1);        % 1 x K
    P_MUSIC = 1 ./ max(denom, eps);

    % 5) 可选平滑与归一化
    if opts.smooth_win > 1
        P_MUSIC = movmean(P_MUSIC, opts.smooth_win);
    end
    if opts.normalize
        P_MUSIC = P_MUSIC ./ max(P_MUSIC);
    end

    % 6) 选峰
    % 将最小峰间距从“度”映射到“索引”
    if numel(theta_scan) > 1
        dth = abs(theta_scan(2) - theta_scan(1));
    else
        dth = 1;
    end
    minDist = max(round(opts.minpkdist_deg / max(dth, eps)), 1);

    if isempty(opts.prominence)
        [~, locs] = findpeaks(P_MUSIC, 'SortStr','descend', ...
                              'MinPeakDistance', minDist);
    else
        [~, locs] = findpeaks(P_MUSIC, 'SortStr','descend', ...
                              'MinPeakDistance', minDist, ...
                              'MinPeakProminence', opts.prominence);
    end
    locs = locs(1:min(nSources, numel(locs)));
    angles_est = theta_scan(locs);
end
