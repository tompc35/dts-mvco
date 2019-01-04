% script T_dens_regression
% ------------------------
% Compute linear regression between temperature and density at site C.

load_dts_isle_data
dens0C_isle = sw_dens0(saltC_isle,wtC_isle);

% regression of density from T
gi = find(isfinite(dens0C_isle));
p = polyfit(wtC_isle(gi),dens0C_isle(gi),1);

% thermal expansion coefficient
alpha = -p(1)/1025; 

plot(wtC_isle,dens0C_isle)