function y = shift_curve(y, off)  %曲线平移函数
% off 可是标量或与 y 同长度的行/列向量
    if isscalar(off)
        y = y + off;
    else
        y = y + reshape(off, size(y));  % 自动匹配形状
    end
end
