clear all

load_dts_isle_data

[x3,y3] = latlon2xy(lat3,lon3,lat3,lon3);
[x4,y4] = latlon2xy(lat4,lon4,lat3,lon3);
[x5,y5] = latlon2xy(lat5,lon5,lat3,lon3);
[xc,yc] = latlon2xy(lat_isle(3),lon_isle(3),lat3,lon3);

dist34 = ((x3-x4)^2 + (y3-y4)^2)^0.5;
dist5c = ((x5-xc)^2 + (y5-yc)^2)^0.5;

% interpolate seagauge temp at C to DTS time base
tsgc = interp1(mdaysg,wtsg(:,3),datetime);

%%

%%% Specify dates %%%

t1 = datenum('13-Jul-2014 0:00');
t2 = datenum('5-Sep-2014 0:00');

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

ti = ti1:ti2;
tim = tim1:tim2;
tis = tis1:tis2;
tid = tis1:tis2;
tivc = tivc1:tivc2;

close all

yl = [17 22];

fs = 8;

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