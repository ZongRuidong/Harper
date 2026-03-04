function [doa_est, spectrum] = wljadmm_doa(Y, theta_scan, A, Rn_est, params)
lambda1 = params.lambda1;
lambdaJ = params.lambdaJ;
epsJ = params.epsJ;
epsL = params.epsL;
r_energy = params.r_energy;  % ✅ 添加这一行
rho = params.rho;
iters = params.iters;
mm_outer = params.mm_outer;
L = params.L;
Rn_est = (Rn_est + Rn_est') / 2;
[Ue, Se] = eig(Rn_est, 'vector');
W = Ue * diag(1 ./ sqrt(max(Se, 1e-9))) * Ue';
Yw = W * Y;
Aw = W * A;


% Step 2: SVD截断
[Uw, Sw, ~] = svd(Yw, 'econ');
s = diag(Sw);
r = find(cumsum(s.^2)/sum(s.^2) >= r_energy, 1);
r = max(1, min(r, size(Uw, 2)));
Ur = Uw(:, 1:r);
Yr = Ur' * Yw;
B = Ur' * Aw;


% Step 3: 列归一化
bn = vecnorm(B, 2, 1) + 1e-12;
Bn = B ./ bn;


% Step 4: 热启动 + MM + ADMM
Z = zeros(size(Bn, 2), size(Yr, 2));
U_admm = Z;
X = Z;
BtB = Bn' * Bn;
BtY = Bn' * Yr;
I = eye(size(BtB));


% 初始解（组Lasso）
for i = 1:iters
X = (BtB + rho * I) \ (BtY + rho * (Z - U_admm));
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
U_admm = U_admm + (X - Z);
end


% MM外层循环
for t = 1:mm_outer
row_norm = sqrt(sum(abs(Z).^2, 2)) + 1e-15;
wJ = 1 ./ (epsJ + row_norm);
wL = 1 ./ (epsL + row_norm);
tau = (lambda1 * wL + lambdaJ * wJ) / rho;


% ADMM再优化
for i = 1:iters
X = (BtB + rho * I) \ (BtY + rho * (Z - U_admm));
X_hat = X + U_admm;
for g = 1:size(X_hat, 1)
row = X_hat(g, :);
norm_row = norm(row, 2);
if norm_row > tau(g)
Z(g, :) = (1 - tau(g)/norm_row) * row;
else
Z(g, :) = 0;
end
end
U_admm = U_admm + (X - Z);
end
end


% Step 5: 恢复尺度 + 计算谱
Xeq = Z ./ bn.';
spectrum = sqrt(sum(abs(Xeq).^2, 2));
spectrum = spectrum(:).';


% Step 6: DOA估计
[~, locs] = findpeaks(spectrum, 'SortStr', 'descend');
doa_est = sort(theta_scan(locs(1:L)));
end
