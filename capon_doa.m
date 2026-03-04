function P_capon = capon_doa(R, A)
    % Capon / MVDR DOA Spectrum Estimation
    % R: 协方差矩阵 MxM
    % A: 导向矩阵 MxK（每列为单位导向矢量）
    
    M = size(R, 1);
    R_inv = pinv(R);  % 可使用 inv 或 pinv
    K = size(A, 2);
    P_capon = zeros(1, K);

    for k = 1:K
        a = A(:, k);
        denom = real(a' * R_inv * a);
        P_capon(k) = 1 / denom;
    end
end
