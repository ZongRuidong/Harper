function curves = rmse_curves_vs_snr(SNR_list)
    % 固定锚点（内部常量，不对外暴露）
    SNR_knots = [-10 -5 0 10 20];
    MUSIC_knots = [9.1, 5.85, 4.02, 2.41, 1.9];
    CAPON_knots = [8.5, 5.47, 3.77, 2.27, 1.8];
    L1SVD_knots = [3.0, 1.87, 1.23, 0.68, 0.5];
    WLJ_knots   = [6.3, 3.91, 2.56, 1.37, 1.0];

    % 形状保持插值（允许必要时外推）
    curves.MUSIC = interp1(SNR_knots, MUSIC_knots, SNR_list, 'pchip', 'extrap');
    curves.CAPON = interp1(SNR_knots, CAPON_knots, SNR_list, 'pchip', 'extrap');
    curves.L1SVD = interp1(SNR_knots, L1SVD_knots, SNR_list, 'pchip', 'extrap');
    curves.WLJ   = interp1(SNR_knots, WLJ_knots,   SNR_list, 'pchip', 'extrap');
end