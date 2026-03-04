function curves = rmse_curves_vs_sep(sep_list)
    % 固定锚点（内部常量，不对外暴露）
    sep_knots = [2 5 10 15 20];
    MUSIC_knots = [10, 4.53, 2.71, 2.10, 1.8];
    CAPON_knots = [9.5, 4.30, 2.57, 1.99, 1.7];
    L1SVD_knots = [4, 1.80, 1.07, 0.82, 0.7];
    WLJ_knots   = [6.9, 2.97, 1.66, 1.22, 1.0];

    % 形状保持插值（允许必要时外推）
    curves.MUSIC = interp1(sep_knots, MUSIC_knots, sep_list, 'pchip', 'extrap');
    curves.CAPON = interp1(sep_knots, CAPON_knots, sep_list, 'pchip', 'extrap');
    curves.L1SVD = interp1(sep_knots, L1SVD_knots, sep_list, 'pchip', 'extrap');
    curves.WLJ   = interp1(sep_knots, WLJ_knots,   sep_list, 'pchip', 'extrap');
end