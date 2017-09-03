function z = zscore_mod(x)

xm = median(x);
mad = median(abs(x-xm));

z = 0.6745*(x-xm)/mad;