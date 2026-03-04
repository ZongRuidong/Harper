function y = adjust_sr_curve(y, cfg, xgrid)
% 与你之前的 adjust_sr_curve 逻辑相同（缩放/平移→向锚点收缩→单调非减→限幅）
    y = y(:).';
    if isfield(cfg,'gain'),  y = y * cfg.gain; end
    if isfield(cfg,'shift'), y = y + cfg.shift; end

    % 逐点锚定并按 mix 融合
    anc = reshape(cfg.anchor, size(y));
    mix = cfg.mix;
    y = mix .* y + (1-mix) .* anc;

    % 限幅与（可选）平滑
    y = min(max(y,0),100);
    if isfield(cfg,'smooth') && cfg.smooth>1
        y = movmean(y, cfg.smooth);
    end

    % 单调非减（保持曲线自然向上）
    if ~isfield(cfg,'enforce') || cfg.enforce
        for i = 2:numel(y)
            y(i) = max(y(i), y(i-1));
        end
    end
    y = reshape(y, size(y));
end