clear all

load_dts_isle_data

[x3,y3] = latlon2xy(lat3,lon3,lat3,lon3);
[x4,y4] = latlon2xy(lat4,lon4,lat3,lon3);
[x5,y5] = latlon2xy(lat5,lon5,lat3,lon3);
[xc,yc] = latlon2xy(lat_isle(3),lon_isle(3),lat3,lon3);

dist34 = ((x3-x4)^2 + (y3-y4)^2)^0.5;
dist5c = ((x5-xc)^2 + (y5-yc)^2)^0.5;

% interpolate seagauge temp at C to DTS time base
nfilt = 3;

% interpolate seagauge temp at C to DTS time base
tsgc = interp1(mdaysg,wtsg(:,3),datetime);
[eventi,event_daten] = get_event_indices_dTdt(boxfilt(tsgc,nfilt),datetime,-0.25);
[eventi3,event_daten3] = get_event_indices_dTdt(boxfilt(tcal3,nfilt),datetime,-0.25);

%% dT/dx and dT/dt

dx = distance(2)-distance(1); % meters
dt = (datetime(2) - datetime(1))*24; % seconds

dTdx = nan(size(tempC));
dTdt = nan(size(tempC));

tempCf = nan(size(tempC));
tempCf(:,2:end-1) = boxfilt(tempC',nfilt)';

dTdx(2:end-1,:) = 0.5*(tempCf(3:end,:)-tempCf(1:end-2,:))/dx;
dTdt(:,2:end-1) = 0.5*(tempCf(:,3:end)-tempCf(:,1:end-2))/dt;

%% dT/dt at moorings
tcf = nan(size(tsgc));
tcf(2:end-1) = boxfilt(tsgc,nfilt);
dTdtc = nan(size(tcf));
dTdtc(2:end-1) = (tcf(3:end)-tcf(1:end-2))/(2*dt);

t3f = nan(size(tcal3));
t3f(2:end-1) = boxfilt(tcal3,nfilt);
dTdt3 = nan(size(t3f));
dTdt3(2:end-1) = (t3f(3:end)-t3f(1:end-2))/(2*dt);

%% Statistics along DTS cable

% Matrix of 1's and 0's indicating presence of front or not
frontpos = zeros(size(tempC));

count_cool = nan(length(distance));
count_warm = nan(length(distance));

dTdt_std = nan(length(distance));

for ci = 1:length(distance)
    Tf_tmp = boxfilt(tempC(ci,:),nfilt);
    dTdt_tmp = dTdt(ci,:);
    
    [eventi_tmp_cool] = get_event_indices_dTdt(Tf_tmp,datetime,-0.25);
    [eventi_tmp_warm] = get_event_indices_dTdt(Tf_tmp,datetime,0.25);
    
    count_cool(ci) = length(eventi_tmp_cool);
    count_warm(ci) = length(eventi_tmp_warm);
    
    frontpos(ci,eventi_tmp_cool) = -1;
    frontpos(ci,eventi_tmp_warm) = 1;
    
    gi = find(isfinite(dTdt_tmp));
    dTdt_std(ci) = std(dTdt_tmp(gi));
end

%% Moving stds

fipr = find(isfinite(prsg(:,3)));
prstd = nan(size(prsg(:,3)));
prstd(fipr) = movingstd(prsg(fipr,3),6*24*2);

fic = find(isfinite(dTdtc));
dTdtc_runstd = nan(size(dTdtc));
dTdtc_runstd(fic) = movingstd(dTdtc(fic),6*24*2);

fi3 = find(isfinite(dTdt3));
dTdt3_runstd = nan(size(dTdt3));
dTdt3_runstd(fi3) = movingstd(dTdt3(fi3),6*24*2);

%%

%%% Specify dates %%%

t1 = datenum('13-Jul-2014 0:00');
t2 = datenum('5-Sep-2014 0:00');

[~,ti1] = min(abs(datetime-t1));
[~,ti2] = min(abs(datetime-t2));
ti = ti1:ti2;

close all

yl = [17 22];

fs = 8;

figure
set(gcf,'paperposition',[0 0 7 4])

subplot(2,1,1)
pcolor(datetime(ti),distance/1000,dTdx(:,ti)), shading flat
xlim([t1 t2])
datetick('x','keeplimits')
set(gca,'tickdir','out')
cbar = colorbar('east');
cpos = get(cbar,'position');
hold on
plot([t1 t2],[distance(zic)/1000 distance(zic)/1000],'k--')
hold on
set(cbar,'position',[cpos(1)+0.1 cpos(2)-0.02 0.5*cpos(3) 0.5*cpos(4)+0.1])
colormap(jet)
ylabel('distance, d [km]','fontsize',fs)
set(gca,'ytick',[1:1:4],'fontsize',fs)
caxis([-0.01,0.01])
colormap(redblue)

subplot(2,1,2)
pcolor(datetime(ti),distance/1000,dTdt(:,ti)), shading flat
xlim([t1 t2])
datetick('x','keeplimits')
set(gca,'tickdir','out')
cbar = colorbar('east');
cpos = get(cbar,'position');
hold on
plot([t1 t2],[distance(zic)/1000 distance(zic)/1000],'k--')
hold on
set(cbar,'position',[cpos(1)+0.1 cpos(2)-0.02 0.5*cpos(3) 0.5*cpos(4)+0.1])
colormap(jet)
ylabel('distance, d [km]','fontsize',fs)
set(gca,'ytick',[1:1:4],'fontsize',fs)
caxis([-1e-3,1e-3])
colormap(redblue)


%%

figure
set(gcf,'paperposition',[0 0 7 4])

subplot(2,1,1)
pcolor(datetime(ti),distance/1000,tempC(:,ti)), shading flat
xlim([t1 t2])
datetick('x','keeplimits')
set(gca,'tickdir','out')
cbar = colorbar('east');
cpos = get(cbar,'position');
hold on
plot([t1 t2],[distance(zic)/1000 distance(zic)/1000],'k--')
hold on
set(cbar,'position',[cpos(1)+0.1 cpos(2)-0.02 0.5*cpos(3) 0.5*cpos(4)+0.1])
colormap(jet)
ylabel('distance, d [km]','fontsize',fs)
set(gca,'ytick',[1:1:4],'fontsize',fs)
title('DTS temperature [deg C]')
xl = xlim;
yld = ylim;
text(xl(1)+0.01*diff(xl),yld(2)-0.073*diff(yld),'a)')

subplot(2,1,2)
thresh = nanmean(tcal3(ti)-tcal4(ti))-2*nanstd(tcal3(ti)-tcal4(ti));
plot(datetime(ti),(tcal3(ti)-tcal4(ti))/dist34,'b-','linewidth',1)
xlim([t1 t2])
hold on
plot(datetime(ti),(tcal5(ti)-tsgc(ti))/dist5c,'r-','linewidth',1)
hold off
legend('R3-R4','R5-C')
datetick('x','keeplimits')
xl = xlim;
yld = ylim;
text(xl(1)+0.01*diff(xl),yld(2)-0.073*diff(yld),'b)')
title('alongshore \DeltaT/\Deltax at bottom')
ylabel('[^oC km^{-1}]','fontsize',fs)
set(gca,'ytick',[-2:1:2],'fontsize',fs)
box on

print -dpng -cmyk -r500 ../figures/fig_event_threshold

%% 
figure

subplot(221)
plot(dTdt_std,distance/1000,'b')
xl = xlim;
hold on
plot(xl,[distance(zic)/1000,distance(zic)/1000],'k--','linewidth',2)
xlabel('std(\partialT/\partialt) [^oC/hr]')
ylabel('distance, d [km]')
xlim(xl);

subplot(222)
hc = plot(count_cool,distance/1000,'b');
hold on
hw = plot(count_warm,distance/1000,'r');
legend([hc(1),hw(1)],'cooling','warming','location','southwest')
xlabel('# of events')
xl = xlim;
plot(xl,[distance(zic)/1000,distance(zic)/1000],'k--','linewidth',2)

subplot(2,2,[3:4])
plot(datetime,dTdt3,'b-')
hold on
plot(datetime,dTdtc,'r-')
datetick('x','keeplimits')
xlim([t1 t2])



figure
subplot(4,1,1)
plot(datetime,dTdt3,'b-')
hold on
plot(datetime,dTdtc,'r-')
yl = ylim;
% plot(datetime(eventi3),yl(1)+0.05*diff(yl),'bx')
% plot(datetime(eventi),yl(1)+0.025*diff(yl),'ro')
xlim([t1 t2])
datetick('x','keeplimits','keepticks')

subplot(4,1,2)
plot(datetime,dTdt3_runstd,'b-')
hold on
plot(datetime,dTdtc_runstd,'r-')
hold off
xlim([t1 t2])
datetick('x','keeplimits','keepticks')

stratC = -(wtC_isle(:,1)-wtC_isle(:,end))/(zsC_isle(1)-zsC_isle(end));
stratH = -(wtH_isle(:,1)-wtH_isle(:,end))/(zsH_isle(1)-zsH_isle(end));

subplot(4,1,3)
plot(mday_isle,stratC,'b-')
hold on
plot(mday_isle,stratH,'r-')
xlim([t1 t2])
datetick('x','keeplimits','keepticks')
ylim([0,0.5])

subplot(4,1,4)
plot(mdaysg,prstd*0.0001,'b-')
xlim([t1 t2])
datetick('x','keeplimits','keepticks')