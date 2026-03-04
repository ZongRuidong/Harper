function y = rescale_powerlaw(x, y1, yN, p)
    z = (x / x(1)).^p;       
    a = (yN - y1) / (z(end) - z(1));
    b = y1 - a * z(1);
    y = a * z + b;
end