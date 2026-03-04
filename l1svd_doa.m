function [doa_est, spectrum] = l1svd_doa(Y, theta_scan, A, params)
% L1-SVD DOA Estimation (不调用外部函数，完整展开)
% Inputs:
% Y - 接收信号 (M x N)
% theta_scan - 扫描角度 (1 x K)
% A - 导向矩阵 (M x K)
% params - 结构体，包含 lambda1, keep_energy, rho, iters, L


lambda1 = params.lambda1;
keep_energy = params.keep_energy;
rho = params.rho;
iters = params.iters;
L = params.L;


% SVD截断
[U, S, ~] = svd(Y, 'econ');
s = diag(S);
r = find(cumsum(s.^2)/sum(s.^2) >= keep_energy, 1);
r = max(1, min(r, size(U, 2)));
Ur = U(:, 1:r);
Yr = Ur' * Y;
B = Ur' * A;


% 列归一化
bn = vecnorm(B, 2, 1) + 1e-12;
Bn = B ./ bn;


% ADMM求解 Z（组Lasso）
Z = zeros(size(Bn, 2), size(Yr, 2));
U_admm = Z;
X = Z;
BtB = Bn' * Bn;
BtY = Bn' * Yr;
I = eye(size(BtB));


for i = 1:iters
% 更新 X
X = (BtB + rho * I) \ (BtY + rho * (Z - U_admm));


% 更新 Z（逐行soft-threshold）
X_hat = X + U_admm;
for g = 1:size(X_hat, 1)
row = X_hat(g, :);
norm_row = norm(row, 2);
if norm_row > lambda1 / rho
Z(g, :) = (1 - lambda1 / (rho * norm_row)) * row;
else
Z(g, :) = 0;
end
end


% 更新 U
U_admm = U_admm + (X - Z);
end


% 还原尺度 & 计算谱
Xeq = Z ./ bn.';
spectrum = sqrt(sum(abs(Xeq).^2, 2));
spectrum = spectrum(:).';


% DOA估计
[~, locs] = findpeaks(spectrum, 'SortStr', 'descend');
doa_est = sort(theta_scan(locs(1:L)));
end