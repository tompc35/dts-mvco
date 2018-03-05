clear all

load_dts_isle_data; 

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

ti = ti1:ti2;
tim = tim1:tim2;
tis = tis1:tis2;
tid = tis1:tis2;

%%
close all

yl = [16.5 22];

set(gcf,'renderer','zbuffer')
set(gcf,'paperposition',[0 0 7 8])

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
ta1 = datenum('15-Jul-2014 7:40');
plot([ta1 ta1],[16.8 17.4],'k')
plot(ta1,17.4,'k^','markerfacecolor','k','markersize',4)
hold off
datetick('x','keeplimits')
%title('H mooring temperature')
ylabel('[deg C]')
ylim(yl)
%grid on
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
ta2 = datenum('15-Jul-2014 20:30');
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
%grid on
pos = get(s2, 'Position') 
posnew = pos; 
posnew(2) = posnew(2) - 0.15; 
set(s2, 'Position', posnew)
set(gca,'color','none')

s1 = subplot(411)
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
ta1 = datenum('15-Jul-2014 7:30');
plot([ta1 ta1],[-0.5 -0.2],'k')
plot(ta1,-0.2,'k^','markerfacecolor','k','markersize',4)

ta2 = datenum('15-Jul-2014 18:30');
plot([ta2 ta2],[3.7 4.0],'k')
plot(ta2,4.0,'k^','markerfacecolor','k','markersize',4)

print -dpng -r600 ../figures/fig_tseries_dts_isle