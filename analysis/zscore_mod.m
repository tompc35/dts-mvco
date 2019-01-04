function z = zscore_mod(x)

% Compute modified z-score based on median.
%
% Reference:
% Iglewicz, B., and D. Hoaglin (1993), How to detect and handle outliers, 
% in ASQC Basic References in Quality Control: Statistical Techniques, 
% vol. 16, edited by E. F. Mykytka, ASQC Quality Press, Milwaukee, WI.

xm = median(x);
mad = median(abs(x-xm));

z = 0.6745*(x-xm)/mad;