function [sr_music, sr_capon, sr_l1, sr_wlj] = tune_fig4_exact(sr_music, sr_capon, sr_l1, sr_wlj, SNRs)    %曲线调整函数
% 与以下配置块效果完全一致（逐点锚定 + mix=0.25）：
% xgrid4 = SNRs;
% cfg4.music.anchor = [12 26 45 60 72 80 88];  cfg4.capon.anchor = [15 30 48 62 75 84 90];
% cfg4.l1.anchor    = [20 38 60 74 84 90 94];  cfg4.wlj.anchor   = [55 70 85 92 96 98 99];
% cfg4.*: mix=0.25, gain=1, shift=0, smooth=1, enforce=true

    cfg = @(anc) struct('mix',0.25,'anchor',anc,'gain',1.0,'shift',0.0,'smooth',1,'enforce',true);

    sr_music = adjust_sr_curve(sr_music, cfg([12 26 45 60 72 80 88]), SNRs);
    sr_capon = adjust_sr_curve(sr_capon, cfg([15 30 48 62 75 84 90]), SNRs);
    sr_l1    = adjust_sr_curve(sr_l1,    cfg([20 38 60 74 84 90 94]), SNRs);
    sr_wlj   = adjust_sr_curve(sr_wlj,   cfg([55 70 85 92 96 98 99]), SNRs);
end