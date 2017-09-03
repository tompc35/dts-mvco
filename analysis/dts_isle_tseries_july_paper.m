clear all

data_dir = '/Users/tomconnolly/work/Data/';

load([data_dir 'ISLE/microcats/microcats.mat'])
mday_isle = jd2matday(jd_isle);

load([data_dir 'ISLE/seagauge/seagauge.mat'])
mdaysg = jd2matday(jdsg');

wtcalfile = [data_dir 'DTS_cal/WaterTempPro/temppro_cals.mat'];
wt = load(wtcalfile);

C = load([data_dir 'ISLE/ADCP/ISLE_C__burstmeans_06-Mar-2015.mat']);
G = load([data_dir 'ISLE/ADCP/ISLE_G__burstmeans_06-Mar-2015.mat']);

Hg = load([data_dir 'ISLE/Gerbi/asitisleburstmeans_02_42_V01.mat']);
Hg.dnum = Hg.dnum_mean-2/24; % data is in Norway time
%%% Find good bins in H velocity %%%
% surface contamination
ht3d = repmat(Hg.height,[length(Hg.dnum_mean) 1 3]);
pr3d = repmat(Hg.pressure,[1 length(Hg.height) 3]);
badi = find(ht3d > pr3d-2);
Hg.ENU = double(Hg.ENU5);
Hg.ENU(badi) = NaN;
clear ht3d pr3d badi
% fouling near end of deployment???
badi = find(Hg.dnum_mean > datenum('Jan-26-2015'));
Hg.ENU(badi,:,:) = NaN;
clear badi

quadpod_dir = [data_dir 'ISLE/quadpod/'];
quadpod_e_mat = [quadpod_dir 'ISLE_E_quadpod_temp_ave_08172016.mat'];
Equad = load(quadpod_e_mat);

dts_dir = [data_dir 'DTS_nc/'];
dts_nc = [dts_dir '/proc/cal_6/DTSasit_proc.nc'];
distance = ncread(dts_nc,'distance');
datetime = ncread(dts_nc,'datetime');
tempC = ncread(dts_nc,'tempC');
[~,zic] = min(abs(distance - 4350)); % corner

cal_nc = [dts_dir 'DTSasit_cal.nc'];
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

% DTS corner index
[~,zic] = min(abs(distance - 4350)); 

t1 = datenum('14-Jul-2014 0:00');
t2 = datenum('18-Jul-2014 12:00');

[~,ti1] = min(abs(datetime-t1));
[~,ti2] = min(abs(datetime-t2));

[~,tim1] = min(abs(mday_isle-t1));
[~,tim2] = min(abs(mday_isle-t2));

[~,tis1] = min(abs(mdaysg-t1));
[~,tis2] = min(abs(mdaysg-t2));

[~,tid1] = min(abs(datetime-t1));
[~,tid2] = min(abs(datetime-t2));

[~,tivc1] = min(abs(C.M.mtime-t1));
[~,tivc2] = min(abs(C.M.mtime-t2));

[~,tivh1] = min(abs(Hg.dnum-t1));
[~,tivh2] = min(abs(Hg.dnum-t2));

ti = ti1:ti2;
tim = tim1:tim2;
tis = tis1:tis2;
tid = tis1:tis2;
tivc = tivc1:tivc2;
tivh = tivh1:tivh2;

%%
close all

yl = [16.5 22];

set(gcf,'renderer','zbuffer')
set(gcf,'paperposition',[0 0 7 8])
%set(gcf,'paperposition',[0 0 8.3 11])

subplot(411)
plot(mday_isle(tim),wtC_isle(tim,:))
xlim([t1 t2])
hold on
plot(mdaysg(tis),wtsg(tis,3),'k-','linewidth',1)
hold off
datetick('x','keeplimits')
title('C mooring temperature')
ylabel('[deg C]')
ylim(yl);
xl = xlim;
text(xl(1)+0.01*diff(xl),yl(2)-0.07*diff(yl),'a)')
set(gca,'ytick',[yl(1):1:yl(2)])
%grid on

subplot(412)
pcolor(datetime(ti),distance/1000,tempC(:,ti)), shading flat
xlim([t1 t2])
datetick('x','keeplimits')
set(gca,'tickdir','out')
cbar = colorbar('east');
cpos = get(cbar,'position');
hold on
plot([t1 t2],[distance(zic)/1000 distance(zic)/1000],'k--')
hold on
yld = ylim;
xl = xlim;
text(xl(1)+0.01*diff(xl),yld(2)-0.07*diff(yld),'b)')
set(cbar,'position',[cpos(1)+0.1 cpos(2)-0.02 0.5*cpos(3) 0.5*cpos(4)+0.1])
colormap(jet)
ylabel('distance, d [km]')
title('DTS cable temperature [deg C]')

subplot(413)
plot(mday_isle(tim),wtH_isle(tim,:))
xlim([t1 t2])
hold on
plot(datetime(ti),tcal3(ti),'k-','linewidth',1)
hold off
datetick('x','keeplimits')
title('H mooring temperature')
ylabel('[deg C]')
ylim(yl)
%grid on
xl = xlim;
text(xl(1)+0.01*diff(xl),yl(2)-0.07*diff(yl),'c)')


subplot(414)
plot(mday_isle(tim),wtE_isle(tim,:))
hold on
plot(Equad.mtime,Equad.aT)
plot(Equad.mtime,Equad.aT(:,1),'k-','linewidth',1)
hold off
xlim([t1 t2])
datetick('x','keeplimits')
title('E mooring temperature')
ylabel('[deg C]')
ylim(yl);
set(gca,'ytick',[yl(1):1:yl(2)])
xl = xlim;
text(xl(1)+0.01*diff(xl),yl(2)-0.07*diff(yl),'d)')

print -dpng -r600 figures_paper/fig_tseries_dts_isle