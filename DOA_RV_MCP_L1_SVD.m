function [DOA_est, P_spatial] = DOA_RV_MCP_L1_SVD(Y, K, theta_grid, d_lambda)
% DOA_RV_MCP_L1_SVD  基于实值化变换与MCP罚函数的L1-SVD算法
%
% 输入:
%   Y          : 接收数据矩阵 [M x T] (复数)
%   K          : 信源数量 (整数)
%   theta_grid : 扫描网格向量 (例如 -90:0.1:90)
%   d_lambda   : 阵元间距与波长比 (通常为 0.5)
%
% 输出:
%   DOA_est    : 估计的角度值 [1 x K]
%   P_spatial  : 归一化空间谱 [1 x N_grid]
%
% 作者: Assistant
% 版本: 1.0 (集成自适应参数与字典归一化)

    %% 1. 基础参数提取
    [M, ~] = size(Y);
    N_theta = length(theta_grid);
    
    %% 2. 实值化预处理 (Data Real-valued Transformation)
    % 2.1 构造稀疏酉矩阵 Q
    Q = loc_construct_Q(M);
    
    % 2.2 变换到波束域
    Y_t = Q' * Y;
    
    % 2.3 实虚分离与扩充 (M x 2T) -> 纯实数
    Y_bar = [real(Y_t), imag(Y_t)];
    
    % 2.4 SVD 降维 (保留信号子空间)
    [U_svd, S_svd, ~] = svd(Y_bar, 'econ');
    % 降维后的观测矩阵 [M x K]
    Y_R = U_svd(:, 1:K) * S_svd(1:K, 1:K);
    
    %% 3. 实数超完备字典构造 (Dictionary Construction)
    % === [关键] 必须使用中心对称索引，否则会有2-3度误差 ===
    p = (-(M-1)/2 : (M-1)/2)'; 
    
    % 生成复数流形
    A_big = exp(1j * 2 * pi * d_lambda * p * sind(theta_grid));
    
    % 变换为实数字典
    B = real(Q' * A_big);
    
    % === [关键] 字典列归一化 (防止虚假峰) ===
    % 计算每一列的模长并归一化
    B_norms = sqrt(sum(B.^2, 1));
    B = bsxfun(@rdivide, B, B_norms);
    
    %% 4. 参数自适应配置
    % 计算最大相关性
    corr_matrix = abs(B' * Y_R);
    max_corr = max(corr_matrix(:));
    
    % Lambda: 稀疏度控制 (经验值: 0.3~0.6 * max_corr)
    lambda = 0.01 * max_corr;
    
    % Gamma: MCP 形状参数
    gamma = 3.7;
    
    % Rho: ADMM 罚参数 (需满足 rho > 1/gamma)
    rho = max(1.5 * lambda, 1.1/gamma + 0.1);
    
    %% 5. 调用核心求解器 (ADMM)
    X_hat = loc_solver_rv_mcp_admm(Y_R, B, lambda, gamma, rho);
    
    %% 6. 谱峰搜索与输出
    % 计算空间谱 (行范数)
    P_raw = sqrt(sum(X_hat.^2, 2))'; % 转为行向量
    
    % 归一化
    P_spatial = P_raw / max(P_raw);
    
    % 简单的谱峰搜索 (Find Peaks)
    [pks, locs] = findpeaks(P_spatial);
    
    % 排序取前 K 个
    if isempty(locs)
        DOA_est = [];
    else
        [~, I] = sort(pks, 'descend');
        if length(I) >= K
            idx_top = locs(I(1:K));
        else
            idx_top = locs(I);
        end
        DOA_est = sort(theta_grid(idx_top));
    end
end

%% ==========================================================
%  以下为本地辅助函数 (Local Functions)
% ==========================================================

function Q = loc_construct_Q(M)
    % 构造稀疏酉矩阵
    eye_k = eye(floor(M/2));      
    J_k = fliplr(eye(floor(M/2))); 
    
    if mod(M, 2) == 0 % 偶数
        Q = (1/sqrt(2)) * [eye_k,       1j*eye_k; 
                           J_k,        -1j*J_k];
    else % 奇数
        row_zeros = zeros(1, floor(M/2));
        col_zeros = zeros(floor(M/2), 1);
        Q = (1/sqrt(2)) * [eye_k,   col_zeros,   1j*eye_k;
                           row_zeros, sqrt(2),     row_zeros;
                           J_k,     col_zeros,  -1j*J_k];
    end
end

function X = loc_solver_rv_mcp_admm(Y_R, B, lambda, gamma, rho)
    % ADMM Solver for RV-MCP
    [~, N_theta] = size(B);
    [~, K] = size(Y_R);
    
    X = zeros(N_theta, K);
    Z = zeros(N_theta, K);
    U = zeros(N_theta, K);
    
    % 预计算矩阵逆 (实数矩阵，速度快)
    % G = (B'B + rho*I)^-1
    G = inv(B' * B + rho * eye(N_theta));
    BtY = B' * Y_R;
    
    max_iter = 150; % 迭代次数
    tol = 1e-5;
    
    for k = 1:max_iter
        % Step 1: Update X
        X = G * (BtY + rho * (Z - U));
        
        % Step 2: Update Z (Firm Thresholding)
        V = X + U;
        Z_old = Z;
        
        % 计算行范数
        v_norms = sqrt(sum(V.^2, 2));
        
        % 区域判断
        idx_zero = v_norms <= lambda/rho;
        idx_pass = v_norms > gamma*lambda;
        idx_scale = (~idx_zero) & (~idx_pass);
        
        % 执行更新
        Z(idx_zero, :) = 0;
        Z(idx_pass, :) = V(idx_pass, :);
        
        if any(idx_scale)
            v_sel = v_norms(idx_scale);
            denom = rho - 1/gamma;
            if denom < 1e-6, denom = 1e-6; end
            scaler = (rho * v_sel - lambda) ./ (v_sel * denom);
            Z(idx_scale, :) = bsxfun(@times, V(idx_scale, :), scaler);
        end
        
        % Step 3: Update U
        U = U + X - Z;
        
        % Check convergence
        if norm(X - Z, 'fro') < tol
            break;
        end
    end
end