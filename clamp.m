function v = clamp(v, lo, hi)  %数值钳制函数
    v = min(max(v,lo),hi);
end