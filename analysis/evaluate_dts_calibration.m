% script calibrate_dts_linearsolve_towercal_threepoints.m
% -------------------------------------------------------
% DTS calibration using full (time-averaged) data set

clear all
run ../data_paths

config_name = 'asit';

dts_nc = [dts_dir 'DTS' config_name '_timeavg.nc'];
cal_nc = [dts_dir 'DTS' config_name '_cal.nc'];
distance = nc_varget(dts_nc,'distance');
datetime = nc_varget(dts_nc,'datetime');

if nc_isvar(dts_nc,'datetime_offset')
    datetime_offset = nc_varget(dts_nc,'datetime_offset');
    datetime = datetime - datetime_offset/24; % adjust all times to GMT
end

wt = load(wtcalfile);

%%

% specify indices of the reference sections
% ~11m offset between distances from splice and distance detected by
% instrument

refzi1 = 15:58;     % cal bath #1 (cooler)
refzi2 = 65:92;     % cal bath #2 (tote)
refzi3 = find(distance >= 132 & distance <= 154);   % anchor / SBE39 0650
refzi4 = find(distance >= 920 & distance <= 930); % WT Pro - 1269446 
refzi5 = find(distance >= 3878 & distance <= 4004); % WT Pro - 1269445
refzi6 = find(distance >= 4572 & distance <= 4582); % WT Pro - 1269447 

datetime = ncread(dts_nc,'datetime');

%%% Load calibration data %%%
tcal1s = ncread(cal_nc,'t_687'); % Seabird
tcal1b = ncread(dts_nc,'tref_2'); % CTEMPS sensor
hi = find(isfinite(tcal1s+tcal1b));
p = polyfit(tcal1b(hi),tcal1s(hi),1);
tcal1_itrp = p(1)*tcal1b + p(2);
tcal1 = tcal1s;
ni = find(isnan(tcal1s));
tcal1(ni) = tcal1_itrp(ni);

tcal2 = ncread(cal_nc,'t_648'); % Seabird
tcal2b = ncread(dts_nc,'tref_1'); % CTEMPS sensor

tcal3 = ncread(cal_nc,'t_650'); % Seabird

% Water Temp Pros (apply calibration)
tcal4p = ncread(cal_nc,'t_446');
j = find(wt.id==1269446);
tcal4=polyval(wt.Tcal(j(end),:),tcal4p);

tcal5p = ncread(cal_nc,'t_445');
j = find(wt.id==1269445);
tcal5=polyval(wt.Tcal(j(end),:),tcal5p);

tcal6p = ncread(cal_nc,'t_447');
j = find(wt.id==1269447);
tcal6=polyval(wt.Tcal(j(end),:),tcal6p);

tint = ncread(dts_nc,'tref_int');

%%% Load DTS data %%%
z1 = mean(ncread(dts_nc,'distance',refzi1(1),length(refzi1)));
stokes1 = ncread(dts_nc,'Stokes',[refzi1(1) 1],[length(refzi1) Inf])';
antistokes1 = ncread(dts_nc,'AntiStokes',[refzi1(1) 1],[length(refzi1) Inf])';
lpr1 = mean(log(stokes1.*antistokes1.^-1),2); % log power ratio, averaged over reference section

z2 = mean(ncread(dts_nc,'distance',refzi2(1),length(refzi2)));
stokes2 = ncread(dts_nc,'Stokes',[refzi2(1) 1],[length(refzi2) Inf])';
antistokes2 = ncread(dts_nc,'AntiStokes',[refzi2(1) 1],[length(refzi2) Inf])';
lpr2 = mean(log(stokes2.*antistokes2.^-1),2); % log power ratio, averaged over reference section

z3 = mean(ncread(dts_nc,'distance',refzi3(1),length(refzi3)));
stokes3 = ncread(dts_nc,'Stokes',[refzi3(1) 1],[length(refzi3) Inf])';
antistokes3 = ncread(dts_nc,'AntiStokes',[refzi3(1) 1],[length(refzi3) Inf])';
lpr3 = mean(log(stokes3.*antistokes3.^-1),2); % log power ratio, averaged over reference section

z4 = mean(ncread(dts_nc,'distance',refzi4(1),length(refzi4)));
stokes4 = ncread(dts_nc,'Stokes',[refzi4(1) 1],[length(refzi4) Inf])';
antistokes4 = ncread(dts_nc,'AntiStokes',[refzi4(1) 1],[length(refzi4) Inf])';
lpr4 = mean(log(stokes4.*antistokes4.^-1),2); % log power ratio, averaged over reference section

z5 = mean(ncread(dts_nc,'distance',refzi5(1),length(refzi5)));
stokes5 = ncread(dts_nc,'Stokes',[refzi5(1) 1],[length(refzi5) Inf])';
antistokes5 = ncread(dts_nc,'AntiStokes',[refzi5(1) 1],[length(refzi5) Inf])';
lpr5 = mean(log(stokes5.*antistokes5.^-1),2); % log power ratio, averaged over reference section

z6 = mean(ncread(dts_nc,'distance',refzi6(1),length(refzi6)));
stokes6 = ncread(dts_nc,'Stokes',[refzi6(1) 1],[length(refzi6) Inf])';
antistokes6 = ncread(dts_nc,'AntiStokes',[refzi6(1) 1],[length(refzi6) Inf])';
lpr6 = mean(log(stokes6.*antistokes6.^-1),2); % log power ratio, averaged over reference section

%%

gamma = nan(size(datetime));
C = nan(size(datetime));
deltaalpha = nan(size(datetime));

disp('calculating calibration parameters')
for ii = 1:length(datetime)
   A = [1  -tcal1(ii)-273.15  (tcal1(ii)+273.15)*z1/1000;
        1  -tcal3(ii)-273.15  (tcal3(ii)+273.15)*z3/1000;
        1  -tcal6(ii)-273.15  (tcal6(ii)+273.15)*z6/1000];  

   b =  [(tcal1(ii)+273.15)*lpr1(ii);
         (tcal3(ii)+273.15)*lpr3(ii);
         (tcal6(ii)+273.15)*lpr6(ii)];
     
   % solve A\b and calculate determinant of A matrix
   if isnan(det(A))
        x = [NaN; NaN; NaN];
        d(ii) = NaN;
   else
        x = A\b;
        d(ii) = cond(A);
   end
   
   gamma(ii) = x(1);
   C(ii) = x(2);
   deltaalpha(ii) = x(3);
   
end

%save dts_cal_params gamma C deltaalpha

% use calibration to calculate temperature at reference sections
t3 = gamma.*(lpr3 + C - deltaalpha*z3/1000).^-1 - 273.15;
t4 = gamma.*(lpr4 + C - deltaalpha*z4/1000).^-1 - 273.15;
t5 = gamma.*(lpr5 + C - deltaalpha*z5/1000).^-1 - 273.15;
t6 = gamma.*(lpr6 + C - deltaalpha*z6/1000).^-1 - 273.15;

t3b = gamma.*(lpr3 + C - nanmean(deltaalpha)*z3/1000).^-1 - 273.15;

%%
gamma_thresh = 300;
gi = find(abs(gamma)>gamma_thresh);
fi = find(abs(gamma)<=gamma_thresh);

t4q = t4;
t4q(fi) = NaN;

%%

figure
subplot(2,2,[1 2])
plot(datetime,[t4],'color',[0.5 0.5 0.5]), datetick('x')
hold on
plot(datetime,[tcal4],'r'), datetick('x')
plot(datetime,t4q,'b')
hold off
ylim([15 21])
title('temperature at R4 [deg C]')
leg = legend('DTS','WTPro','DTS (QC)');
set(leg,'box','off')
set(gca,'fontsize',14)
xlim([datenum('July 10 2014') datenum('Oct 20 2014')])
ylim([15 21.5])

subplot(2,2,3)
h =patch([-300 -300 300 300],[0 12 12 0],[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'linestyle','none','edgecolor',[0.5 0.5 0.5]);
hold on
plot(gamma,abs(t4-tcal4),'.','color','b','markersize',8)
hold off
ylim([0 12])
xlim([-4000 4000])
ylabel('T error [deg C]')
xlabel('\gamma [deg K]')
set(gca,'fontsize',14)
box on

subplot(2,2,4)
plot(t4(gi),tcal4(gi),'b.')
xlabel('cable T [deg C]')
ylabel('sensor T [deg C]')
rmse = sqrt(nanmean((t4(gi)-tcal4(gi)).^2));
bias = nanmean(t4(gi)-tcal4(gi));
set(gca,'fontsize',14)
ylim([15 21])
xlim([15 21])
axis square
hold on
plot([15 21],[15 21],'--','color',[0.5 0.5 0.5])
set(gca,'xtick',[15:21])
hold off
yl = ylim;
xl = xlim;
text(xl(1)+diff(xl)*.05,yl(2)-diff(yl)*.05,['RMSE = ' num2str(rmse,2) '^oC'])
text(xl(1)+diff(xl)*.05,yl(2)-diff(yl)*.15,['bias = ' num2str(bias,2) '^oC'])

set(gcf, 'PaperSize', [7 5.5]); %Keep the paper size [width height] 
set(gcf, 'PaperPosition', [0 0 7 5.5]); %

print -dpdf ../figures/fig_DTS_error_summary