function y = lin_endpoints(x, y1, yN)  %线性端点插值函数
    n = numel(x);
    if n == 1
        y = y1;
    else
        t = linspace(0, 1, n);
        y = y1 + (yN - y1) * t;
    end
end