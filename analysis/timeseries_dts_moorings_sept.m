% script timeseries_dts_moorings_sept.m
% -------------------------------------
% Plot time series of DTS and mooring temperatures - 4.5 day subset in
% September

clear all

load_dts_isle_data; 

t1 = datenum('1-Sep-2014 9:00');
t2 = t1+4.5;

[~,ti1] = min(abs(datetime-t1));
[~,ti2] = min(abs(datetime-t2));

[~,tim1] = min(abs(mday_isle-t1));
[~,tim2] = min(abs(mday_isle-t2));

[~,tis1] = min(abs(mdaysg-t1));
[~,tis2] = min(abs(mdaysg-t2));

[~,tid1] = min(abs(datetime-t1));
[~,tid2] = min(abs(datetime-t2));

ti = ti1:ti2;
tim = tim1:tim2;
tis = tis1:tis2;
tid = tis1:tis2;

%%
close all

yl = [16.5 22];

set(gcf,'renderer','zbuffer')
set(gcf,'paperposition',[0 0 7 8])
%set(gcf,'paperposition',[0 0 8.3 11])

subplot(414)
plot(mday_isle(tim),wtE_isle(tim,:))
hold on
plot(Equad.mtime,Equad.aT)
plot(Equad.mtime,Equad.aT(:,1),'k-','linewidth',1)
hold off
xlim([t1 t2])
datetick('x','keeplimits')
ylabel('[deg C]')
ylim(yl);
set(gca,'ytick',[yl(1):1:yl(2)])
xl = xlim;
text(xl(1)+0.01*diff(xl),yl(2)-0.07*diff(yl),'    Site E')
box off
set(gca,'color','none')

s3 = subplot(413);
plot(mday_isle(tim),wtH_isle(tim,:))
xlim([t1 t2])
hold on
plot(datetime(ti),tcal3(ti),'k-','linewidth',1)
ta1 = datenum('2-Sep-2014 11:30');
plot([ta1 ta1],[17.2 17.8],'k')
plot(ta1,17.8,'k^','markerfacecolor','k','markersize',4)
hold off
datetick('x','keeplimits')
ylabel('[deg C]')
ylim(yl)
xl = xlim;
text(xl(1)+0.01*diff(xl),yl(2)-0.07*diff(yl),'    Site H')
box off
set(gca,'xticklabel',{[]}) 
set(gca,'xaxisloc','top')
set(gca,'yaxisloc','right')
set(gca,'xcolor','none')
pos = get(s3, 'Position') 
posnew = pos; 
posnew(2) = posnew(2) - 0.09; 
set(s3, 'Position', posnew)
set(gca,'ytick',15.5:1:30)
set(gca,'color','none')

s2 = subplot(412);
plot(mday_isle(tim),wtC_isle(tim,:))
xlim([t1 t2])
hold on
plot(mdaysg(tis),wtsg(tis,3),'k-','linewidth',1)
ta2 = datenum('3-Sep-2014 03:30');
plot([ta2 ta2],[17.6 18.2],'k')
plot(ta2,18.2,'k^','markerfacecolor','k','markersize',4)
hold off
datetick('x','keeplimits')
title('mooring temperature [deg C]')
ylabel('[deg C]')
ylim(yl);
xl = xlim;
text(xl(1)+0.01*diff(xl),yl(2)+0.11*diff(yl),'b)')
text(xl(1)+0.01*diff(xl),yl(2)-0.07*diff(yl),'    Site C')
set(gca,'ytick',[yl(1):1:yl(2)])
set(gca,'xticklabel',{[]}) 
set(gca,'xaxisloc','top')
set(gca,'xcolor','none')
box off
pos = get(s2, 'Position') 
posnew = pos; 
posnew(2) = posnew(2) - 0.15; 
set(s2, 'Position', posnew)
set(gca,'color','none')

ta1 = datenum('2-Sep-2014 10:30');
ta2 = datenum('3-Sep-2014 01:30');
dt = 1/24/60;
tt = ta1:dt:ta2;
uc = C.M.evm(3,:); 
uh = H.uu(3,:); 
uth = interp1(H.ttime,-uh,tt);
utc = interp1(C.M.mtime,-uc,tt);
utavg = 0.5*(interp1(H.ttime,-uh,tt) + interp1(C.M.mtime,-uc,tt));

s1 = subplot(411)
pcolor(datetime(ti),distance/1000,tempC(:,ti)), shading flat
xlim([t1 t2])
datetick('x','keeplimits')
set(gca,'tickdir','out')
cbar = colorbar('east');
cpos = get(cbar,'position');
hold on
plot(tt,cumtrapz(utc*86400*dt)/1000+distance(1)/1000,'--','color','k','linewidth',2)
text(tt(end),1.82,'C velocity','color','k')
plot(tt(1:450),cumtrapz(uth(1:450)*86400*dt)/1000+distance(1)/1000,'-','color','k','linewidth',2)
text(tt(450),0.23,'H velocity','color','k')
plot([t1 t2],[distance(zic)/1000 distance(zic)/1000],'k--')
hold on
yld = ylim;
xl = xlim;
text(xl(1)+0.01*diff(xl),yld(2)-0.07*diff(yld),'a)')
text(xl(2)+0.01*diff(xl),distance(zic)/1000,'site C')
text(xl(2)+0.01*diff(xl),0,'site H')
set(cbar,'position',[cpos(1)+0.1 cpos(2)-0.09 0.5*cpos(3) 0.5*cpos(4)+0.1])
colormap(jet)
ylabel('distance, d [km]')
title('DTS cable temperature [deg C]')
p = get(s1, 'position');
offset = 0.11;
set(s1, 'position', [p(1), p(2)-offset, p(3), p(4)+offset]);

ax = gca;
ax.Clipping = 'off';
plot([ta1 ta1],[-0.5 -0.2],'k')
plot(ta1,-0.2,'k^','markerfacecolor','k','markersize',4)

plot([ta2 ta2],[3.7 4.0],'k')
plot(ta2,4.0,'k^','markerfacecolor','k','markersize',4)

print -djpeg -r600 ../figures/fig_tseries_dts_isle_sept.jpg