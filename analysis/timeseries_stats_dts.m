clear all

do_corr = 0;

load_dts_isle_data
detide_c_components

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

tsgc_isle = interp1(mdaysg,wtsg(:,3),mday_isle);
t3_isle = interp1(datetime,tcal3,mday_isle);

%%% Specify dates %%%

t1 = datenum('13-Jul-2014 0:00');
t2 = datenum('5-Sep-2014 0:00');

[~,ti1] = min(abs(datetime-t1));
[~,ti2] = min(abs(datetime-t2));
ti = ti1:ti2;

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

t4f = nan(size(tcal4));
t4f(2:end-1) = boxfilt(tcal4,nfilt);
dTdt4 = nan(size(t4f));
dTdt4(2:end-1) = (t4f(3:end)-t4f(1:end-2))/(2*dt);

t5f = nan(size(tcal5));
t5f(2:end-1) = boxfilt(tcal5,nfilt);
dTdt5 = nan(size(t5f));
dTdt5(2:end-1) = (t5f(3:end)-t5f(1:end-2))/(2*dt);

t6f = nan(size(tcal6));
t6f(2:end-1) = boxfilt(tcal6,nfilt);
dTdt6 = nan(size(t6f));
dTdt6(2:end-1) = (t6f(3:end)-t6f(1:end-2))/(2*dt);

%% Statistics along DTS cable

% Matrix of 1's and 0's indicating presence of front or not
frontpos = zeros(size(tempC));

r3 = nan(size(distance));
lag3 = nan(size(distance));

rc = nan(size(distance));
lagc = nan(size(distance));

r4 = nan(size(distance));
lag4 = nan(size(distance));

r5 = nan(size(distance));
lag5 = nan(size(distance));

r6 = nan(size(distance));
lag6 = nan(size(distance));

count_cool = nan(size(distance));
count_warm = nan(size(distance));

dTdt_std = nan(size(distance));

disp('cable event and cross-correlation loop')
for ci = 1:length(distance)
    if rem(ci,10) == 0
       disp([num2str(ci) '/' num2str(length(distance))]) 
    end
    
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
    
    if do_corr
        [rhoxy,alln,dof,rhoxx,rhoyy,rhoyx,inti] = xcov_dof(dTdt3(ti),dTdt(ci,ti)');
        [rmax,ri] = max(rhoxy);
        r3(ci) = rmax;
        lag3(ci) = alln(ri);
        
        [rhoxy,alln,dof,rhoxx,rhoyy,rhoyx,inti] = xcov_dof(dTdt4(ti),dTdt(ci,ti)');
        [rmax,ri] = max(rhoxy);
        r4(ci) = rmax;
        lag4(ci) = alln(ri);
        
        [rhoxy,alln,dof,rhoxx,rhoyy,rhoyx,inti] = xcov_dof(dTdt5(ti),dTdt(ci,ti)');
        [rmax,ri] = max(rhoxy);
        r5(ci) = rmax;
        lag5(ci) = alln(ri);
        
        [rhoxy,alln,dof,rhoxx,rhoyy,rhoyx,inti] = xcov_dof(dTdt6(ti),dTdt(ci,ti)');
        [rmax,ri] = max(rhoxy);
        r6(ci) = rmax;
        lag6(ci) = alln(ri);

        [rhoxy,alln,dof,rhoxx,rhoyy,rhoyx,inti] = xcov_dof(dTdtc(ti),dTdt(ci,ti)');
        [rmax,ri] = max(rhoxy);
        rc(ci) = rmax;
        lagc(ci) = alln(ri);

    end
end

%% Moving stds

navg = 6*24*2;

fipr = find(isfinite(prsg(:,3)));
prstd = nan(size(prsg(:,3)));
prstd(fipr) = movingstd(prsg(fipr,3),navg);

ficm = find(isfinite(C.M.evm(1,:)));
ucstd = nan(size(C.M.evm(1,:)));
ucstd(ficm) = movingstd(C.M.evm(1,ficm),navg);

ficmt = find(isfinite(C.utide(1,:)));
uctstd = nan(size(C.utide(1,:)));
uctstd(ficmt) = movingstd(C.utide(1,ficmt),navg);

ficmh = find(isfinite(H.uu(2,:)));
uhstd = nan(size(H.uu(2,:)));
uhstd(ficmh) = movingstd(H.uu(2,ficmh),navg);

fic = find(isfinite(dTdtc));
dTdtc_runstd = nan(size(dTdtc));
dTdtc_runstd(fic) = movingstd(dTdtc(fic),navg);

fi3 = find(isfinite(dTdt3));
dTdt3_runstd = nan(size(dTdt3));
dTdt3_runstd(fi3) = movingstd(dTdt3(fi3),navg);

%%% Stratification

% density
stratC = (9.8/1025)*(sigC_isle(:,1)-sigC_isle(:,end))/(zsC_isle(1)-zsC_isle(end));
stratH = (9.8/1025)*(sigH_isle(:,1)-sigH_isle(:,end))/(zsH_isle(1)-zsC_isle(end));

stratC_interp = interp1(mday_isle,stratC,datetime);
stratH_interp = interp1(mday_isle,stratH,datetime);

stratH_filt = pl64tc(stratH_interp,6*33);
stratC_filt = pl64tc(stratC_interp,6*33);

stratH_boxfilt = nan(size(stratH_interp));
stratC_boxfilt = nan(size(stratH_interp));
stratH_boxfilt(navg/2:end-navg/2) = boxfilt(stratH_interp,navg);
stratC_boxfilt(navg/2:end-navg/2) = boxfilt(stratC_interp,navg);

% temperature - full water column
wtC_bot = interp1(mdaysg,wtsg(:,3),mday_isle);
wtH_bot = interp1(datetime,tcal3,mday_isle);

stratC_temp = -(wtC_isle(:,1)-wtC_bot(:,end))/(zsC_isle(1)-15);
stratH_temp = -(wtH_isle(:,1)-wtH_bot(:,end))/(zsH_isle(1)-16);

stratC_temp_interp = interp1(mday_isle,stratC_temp,datetime);
stratH_temp_interp = interp1(mday_isle,stratH_temp,datetime);

stratH_temp_filt = pl64tc(stratH_temp_interp,6*33);
stratC_temp_filt = pl64tc(stratC_temp_interp,6*33);

stratH_temp_boxfilt = nan(size(stratH_temp_interp));
stratC_temp_boxfilt = nan(size(stratH_temp_interp));
stratH_temp_boxfilt(144:end-144) = boxfilt(stratH_temp_interp,6*24*2);
stratC_temp_boxfilt(144:end-144) = boxfilt(stratC_temp_interp,6*24*2);

gfi = find(isfinite(stratH_temp_boxfilt(ti)+dTdt3_runstd(ti)));

% correlation
rmat = corrcoef(stratH_temp_boxfilt(ti),dTdt3_runstd(ti));
[rhoxy,alln,dof,rhoxx,rhoyy,rhoyx,inti] = xcov_dof(stratH_temp_boxfilt(ti),dTdt3_runstd(ti));
ro = rmat(2);
p = rsig(ro,dof);

disp('correlation - std(dT/dt) vs. dT/dz')
disp('----------------------------------')
disp(['r = ' num2str(ro)])
disp(['dof = ' num2str(dof)])
disp(['p = ' num2str(p)])

%%

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

%%

figure
plot(stratH_temp_boxfilt(ti(gfi)),dTdt3_runstd(ti(gfi)),'.')

figure
set(gcf,'papersize',[6 7])
set(gcf,'paperposition',[0 -0.4 6 7.6])

subplot(4,1,1)
plot(datetime,dTdt3,'b-')
hold on
plot(datetime,dTdtc,'r-')
yl = ylim;
xlim([t1 t2])
datetick('x','keeplimits','keepticks')
ylabel('\it\fontname{Times}\partialT/\partialt \rm\fontname{Helvetica}[^oC/hr]')
leg = legend('H','C','location','southeast');
set(leg,'box','off')
title('\rm time derivative of bottom temperature, \it\fontname{Times}\partialT/\partialt')
set(gca,'xticklabel',[])
ylim([-4 2])
xl = xlim;
yl = ylim;
text(xl(1)+0.01*diff(xl),yl(2)-0.08*diff(yl),'a)')

subplot(4,1,2)
plot(datetime,dTdt3_runstd,'b-')
hold on
plot(datetime,dTdtc_runstd,'r-')
hold off
ylim([0,0.5])
xlim([t1 t2])
datetick('x','keeplimits','keepticks')
ylabel(['std(\it\fontname{Times}\partialT/\partialt\rm\fontname{Helvetica}) [^oC/hr]'])
title('\rm 2-day running standard deviation of \it\fontname{Times}\partialT/\partialt')
set(gca,'xticklabel',[])
leg = legend('H','C','location','northeast');
set(leg,'box','off')
xl = xlim;
yl = ylim;
text(xl(1)+0.01*diff(xl),yl(2)-0.08*diff(yl),'b)')

subplot(4,1,3)
plot(datetime,stratH_temp_boxfilt,'b-')
hold on
plot(datetime,stratC_temp_boxfilt,'r-')
xlim([t1 t2])
datetick('x','keeplimits','keepticks')
ylim([0,0.25])
set(gca,'xticklabel',[])
title('\rm 2-day running mean of vertical temperature stratification, \it\fontname{Times}\partialT/\partialz')
ylabel('\it\fontname{Times}\partialT/\partialz\rm\fontname{Helvetica} [^oC/m]')
leg = legend('H','C','location','northeast');
set(leg,'box','off')
xl = xlim;
yl = ylim;
text(xl(1)+0.01*diff(xl),yl(2)-0.08*diff(yl),'c)')

subplot(4,1,4)
plot(C.M.mtime,ucstd,'r-')
hold on
plot(C.M.mtime,uctstd,'r--')
xlim([t1 t2])
datetick('x','keeplimits','keepticks')
title('\rm 2-day running standard deviation of eastward velocity (site C)')
ylabel('std(\it\fontname{Times}u\rm\fontname{Helvetica}) [m/s]')
leg = legend('full','M_2+S_2+N_2','location','northeast');
set(leg,'box','off');
ylim([0.07,0.11])
xl = xlim;
yl = ylim;
text(xl(1)+0.01*diff(xl),yl(2)-0.08*diff(yl),'d)')
set(gcf,'renderer','painters')

print -dpdf ../figures/fig_tseries_dTdt_strat_u

%%
figure
set(gcf,'papersize',[6 5])
set(gcf,'paperposition',[-0.3 0 6.5 5])

subplot(221)
plot(dTdt_std,distance/1000,'k')
xl = xlim;
hold on
plot(xl,[distance(zic)/1000,distance(zic)/1000],'k--','linewidth',2)
xlabel('std(\it\fontname{Times}\partialT/\partialt\rm\fontname{Helvetica}) [^oC/hr]')
ylabel('distance, d [km]')
xlim(xl);
yl = ylim;
text(xl(1)+0.01*diff(xl),yl(2)-0.06*diff(yl),'a)')

subplot(222)
hc = plot(count_cool,distance/1000,'b');
hold on
hw = plot(count_warm,distance/1000,'r');
leg = legend([hc(1),hw(1)],'cooling','warming','location','southwest');
set(leg,'box','off')
xlabel('# of events')
xl = xlim;
plot(xl,[distance(zic)/1000,distance(zic)/1000],'k--','linewidth',2)
text(xl(1)+0.01*diff(xl),yl(2)-0.06*diff(yl),'b)')

subplot(223)
[~,zi3] = min(abs(distance-d3));
[~,zi4] = min(abs(distance-d4));
[~,zi5] = min(abs(distance-d5));
[~,zi6] = min(abs(distance-d6));
h = plot(r3,distance/1000); c3 = get(h,'color');
hold on
h = plot(r4,distance/1000); c4 = get(h,'color');
h = plot(rc,distance/1000); cc = get(h,'color');
h = plot(r5,distance/1000); c5 = get(h,'color');
h = plot(r6,distance/1000); c6 = get(h,'color');
plot(r3(zi3),d3/1000,'o','color',c3,'markerfacecolor',c3)
plot(r4(zi4),d4/1000,'o','color',c4,'markerfacecolor',c4)
plot(rc(zic),distance(zic)/1000,'o','color',cc,'markerfacecolor',cc)
plot(r5(zi5),d5/1000,'o','color',c5,'markerfacecolor',c5)
plot(r6(zi6),d6/1000,'o','color',c6,'markerfacecolor',c6)
xl = xlim;
plot(xl,[distance(zic)/1000,distance(zic)/1000],'k--','linewidth',2)
xlabel('correlation coefficient, r')
ylabel('distance, d [km]')
leg = legend('r3','r4','C','r5','r6','location','east');
set(leg,'box','off')
text(xl(1)+0.01*diff(xl),yl(2)-0.06*diff(yl),'c)')

subplot(224)
rcrit = 0.24; % only show lags for correlations above this threshold
plot(lag3(find(r3>rcrit))/6,distance(find(r3>rcrit))/1000)
hold on
plot(lag4(find(r4>rcrit))/6,distance(find(r4>rcrit))/1000)
plot(lagc(find(rc>rcrit))/6,distance(find(rc>rcrit))/1000)
plot(lag5(find(r5>rcrit))/6,distance(find(r5>rcrit))/1000)
plot(lag6(find(r6>0.4))/6,distance(find(r6>0.4))/1000)
plot(lag3(zi3)/6,d3/1000,'o','color',c3,'markerfacecolor',c3)
plot(lag4(zi4)/6,d4/1000,'o','color',c4,'markerfacecolor',c4)
plot(lagc(zic)/6,distance(zic)/1000,'o','color',cc,'markerfacecolor',cc)
plot(lag5(zi5)/6,d5/1000,'o','color',c5,'markerfacecolor',c5)
plot(lag6(zi6)/6,d6/1000,'o','color',c6,'markerfacecolor',c6)
xlim([-1.5,1.5])
xl = xlim;
plot(xl,[distance(zic)/1000,distance(zic)/1000],'k--','linewidth',2)
xlabel('lag [hours]')
text(xl(1)+0.01*diff(xl),yl(2)-0.06*diff(yl),'d)')

print -dpdf ../figures/fig_cable_stats

%%
dcorr3 = max(abs(distance(zi3) - distance(find(r3*exp(1)>r3(zi3)))));
dcorr4 = max(abs(distance(zi4) - distance(find(r4*exp(1)>r4(zi4)))));
dcorr5 = max(abs(distance(zi5) - distance(find(r5*exp(1)>r5(zi5)))));
dcorr6 = max(abs(distance(zi6) - distance(find(r6*exp(1)>r6(zi6)))));
dcorrc = max(abs(distance(zic) - distance(find(rc*exp(1)>rc(zic)))));

disp('decorrelation scales')
disp('--------------------')
disp(dcorr3)
disp(dcorr4)
disp(dcorr5)
disp(dcorr6)
disp(dcorrc)