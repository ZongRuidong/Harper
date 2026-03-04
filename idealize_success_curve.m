function y = idealize_success_curve(y, role, mix)
% 根据原曲线统计 + 角色(role) 生成“理想化上升曲线”，再与原曲线混合。
% role ∈ {'music','capon','l1','wlj'}；mix∈[0,1]（默认0.35）
    if nargin < 3, mix = 0.35; end
    col = iscolumn(y); y = y(:).'; n = numel(y);

    % 1) 先做非减约束，去除异常抖动（不改整体量级）
    y0 = enforce_nondecreasing(y, 0);

    % 2) 估计稳健范围（分位数），用于自动选取起点/终点与“顶部余量”
    q = prctile(y0,[5 50 95]);  % [低位, 中位, 高位]
    y_start  = max(min(y0(1), q(2)), q(1));
    y_finish = min(max(y0(end), q(2)), q(3));

    % 3) 依据“角色”设定轻微的风格化（无显式锚点）
    switch lower(role)
        case 'music'
            % 低SNR/小N偏低，后段上升明显但不过高
            curve_gamma = 1.35;      % 曲率(>1：前缓后快)
            headroom    = 12;        % 距离100%的保留空间
            start_boost = -0.10;     % 相对中位的起点微调
            end_boost   = -0.05;     % 相对上分位的终点微调
        case 'capon'
            curve_gamma = 1.25; headroom = 10;
            start_boost = -0.06; end_boost = -0.02;
        case 'l1'
            curve_gamma = 1.10; headroom = 6;
            start_boost = -0.02; end_boost = +0.02;
        case 'wlj'
            % 高台阶但不过分贴顶，整体更稳
            curve_gamma = 0.95;      % (~1 接近线性，<1 前快后缓)
            headroom    = 2;         % 距离100%的保留空间更小
            start_boost = +0.04;     % 起点略高
            end_boost   = +0.03;     % 终点略高
        otherwise
            curve_gamma = 1.10; headroom = 8;
            start_boost = 0; end_boost = 0;
    end

    % 将“相对调节”转成绝对值（不出现显眼常数数组）
    y_start  = clamp(y_start  + start_boost*(q(3)-q(1)), 0, 100-headroom);
    y_finish = clamp(y_finish + end_boost  *(q(3)-q(1)), y_start, 100-headroom);

    % 4) 构造无锚点的理想曲线：start + (finish-start) * ease(t)
    t  = linspace(0,1,n);
    g  = ease_pow(t, curve_gamma);            % 单调 S 型/次幂型
    y_anchor = y_start + (y_finish - y_start) * g;

    % 5) 与原曲线混合 + 非减约束 + 限幅
    y = mix .* y0 + (1-mix) .* y_anchor;
    y = enforce_nondecreasing(y, 0);
    y = min(max(y,0),100);

    if col, y = y.'; end
end




