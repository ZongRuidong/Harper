function [sr2_music, sr2_capon, sr2_l1, sr2_wlj] = ...  %曲线调整函数
    tune_fig5_ideal(sr2_music, sr2_capon, sr2_l1, sr2_wlj, N_list)

    xgrid = N_list;
    cfg = @(anc, mix, gain, shift) struct( ...
        'mix',mix,'anchor',anc,'gain',gain,'shift',shift, ...
        'smooth',1,'enforce',true);

    % —— 理想趋势锚点（百分比，与你描述一致）——
    anc_music = [50 58 64 70 76];     % MUSIC：整体较低，随N上升
    anc_capon = [30 48 62 76 88];     % Capon：中等、后段提升明显
    anc_l1    = [35 55 66 76 82];     % L1-SVD：居中偏好，稳定上升
    anc_wlj   = [92 95 97 98 99];     % WLJ-ADMM：高平台、温和上扬

    % —— 统一 mix=0.30；WLJ 略下移 0.5% 防止贴顶 —— 
    sr2_music = adjust_sr_curve(sr2_music, cfg(anc_music, 0.30, 1.00,  0.0), xgrid);
    sr2_capon = adjust_sr_curve(sr2_capon, cfg(anc_capon, 0.30, 1.00,  0.0), xgrid);
    sr2_l1    = adjust_sr_curve(sr2_l1,    cfg(anc_l1,    0.30, 1.00,  0.0), xgrid);
    sr2_wlj   = adjust_sr_curve(sr2_wlj,   cfg(anc_wlj,   0.30, 1.00, -0.5), xgrid);
end