clear all

plot_event = 1;

load_dts_isle_data
nfilt = 3;

detide_c;

% interpolate seagauge temp at C to DTS time base
tsgc = interp1(mdaysg,wtsg(:,3),datetime);
[eventi,event_daten] = get_event_indices_dTdt(boxfilt(tsgc,nfilt),datetime,-0.25);

%%

% 5 - calibration point 5
% c - corner of DTS cable
% e - end of DTS cable

ti = 1:length(eventi);

tc = interp1(mday_isle,wtC_isle,datetime);

close all
clear c_*
clear cpx cpy cpax cpay
clear *_events

% indices of events to plot (three)
pii = [1,5,15];

ii = 1;
for kk = 1:length(ti)
    
    jj = ti(kk);
    t1 = event_daten(jj)-0.25;
    t2 = event_daten(jj)+0.25;
    di = find(datetime>=t1 & datetime<=t2);
    Ci = find(C.M.mtime>=t1 & C.M.mtime<t2);
    tim = find(mday_isle>=t1 & mday_isle<=t2);
    
    if ii == 1
        datetime_events = datetime(di)-event_daten(jj);
        mday_isle_events = mday_isle(tim)-event_daten(jj);
        mitime_events = C.M.mtime(Ci)-event_daten(jj);
    end
    
    % use only events where data is good (check corner)
    if ~isnan(tempC(zic,di))
    
        for zii = 1:40
            zjj = zii;
            
            % modify points which 3 points are being used in calculation
            zicr = zic; %  corner
            zi5 = zic-zii;
            zie = zic+zjj;

            latc = lat_dts(zicr);
            lonc = lon_dts(zicr);

            late = lat_dts(zie);
            lone = lon_dts(zie);
            
            lat5 = lat_dts(zi5);
            lon5 = lon_dts(zi5);

            tcrn = tempC(zicr,:);
            tend = tempC(zie,:);
            t5 = tempC(zi5,:);   
            
            tcrn = boxfilt(tcrn,nfilt);
            tend = boxfilt(tend,nfilt);
            t5 = boxfilt(t5,nfilt);
            
            dt = (datetime(2)-datetime(1))*86400;
            
            [cp,phixy,ic,ie,i5] = phase_velocity(lonc,latc,tcrn(di),lone,late,tend(di),lon5,lat5,t5(di),dt);
            
            % find start and end time of event passage
            %tstart = datetime(di(min([i5 ic ie])));
            %tf = datetime(di(max([i5 ic ie])));
            tstart = datetime(di(ic)-7);
            tf = datetime(di(ic)-5);            

            % mean tidal velocity in direction of phase propagation during
            % passage of event
            cmi = find(C.M.mtime>=tstart & C.M.mtime<=tf);
            %cmi = max(cmi);
            if ~isempty(cmi)
                cmih = [(cmi(1)-2):(cmi(1)-1)];
                %cmih = [(cmi(end)+3):(cmi(end)+4)];
                %cmih = [(cmi(end)+1):(cmi(end)+2)];
                %cmih = round(mean(cmi));
            else
                cmih = cmi;
            end

            zoff = 0.75; % difference between bottom and ADCP height

            uhm_tide = mean(depthavg(C.M.evm(1:17,cmi),15-(C.M.z(1:17)+zoff),15));
            vhm_tide = mean(depthavg(C.M.nvm(1:17,cmi),15-(C.M.z(1:17)+zoff),15));
            
            %uhm_tide = mean(C.M.evm(1,cmi));
            %vhm_tide = mean(C.M.nvm(1,cmi));
            
            %uhm_tide = mean(depthavg(C.utide(1:17,cmih),15-(C.M.z(1:17)+zoff),15),2);
            %vhm_tide = mean(depthavg(C.vtide(1:17,cmih),15-(C.M.z(1:17)+zoff),15),2);

            whm_tide = uhm_tide+i*vhm_tide;

            whmc_tide = whm_tide*exp(-i*phixy);
            uhmc_tide = real(whmc_tide);   

            wc = C.M.evm + i*C.M.nvm;
            %wcc = wc*exp(-i*phixy);

            % remove background current
            cpa = cp-mean(uhmc_tide);

            dx_d(zii) = (distance(end)-distance(end-1))*(zii);
            dy_d(zii) = (distance(end)-distance(end-1))*(zjj);

            ur_d(zii) = mean(uhmc_tide);
            u_d(zii) = mean(real(whm_tide));
            v_d(zii) = mean(imag(whm_tide));

            cp_d(zii) = cp;
            cpa_d(zii) = cpa;
            phixy_d(zii) = phixy;

            cpx_d(zii) = cp*cos(phixy);
            cpy_d(zii) = cp*sin(phixy);
            cpax_d(zii) = cpa*cos(phixy);
            cpay_d(zii) = cpa*sin(phixy);
            
            whm_tide_d(zii) = whm_tide;
            uhmc_tide_d(zii) = uhmc_tide;
        end
        
        % choose the range of distances to analyze
        m = 1;
        n = 15;
        
        cpx_dv = vecshape(cpx_d(m:n));
        cpy_dv = vecshape(cpy_d(m:n));
        
        cpax_dv = vecshape(cpax_d(m:n));
        cpay_dv = vecshape(cpay_d(m:n));
        
        fi = find(isfinite(cpx_dv+cpy_dv));
        cpx_dv = cpx_dv(fi);
        cpy_dv = cpy_dv(fi);
        
        fi = find(isfinite(cpax_dv+cpay_dv));
        cpax_dv = cpax_dv(fi);
        cpay_dv = cpay_dv(fi);
        
        % Remover outliers using modified z-score of  Iglewicz and Hoaglin (1993)
        % http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
        gi = find(abs(zscore_mod(cpx_dv)) <= 3.5 & ...
                  abs(zscore_mod(cpy_dv)) <= 3.5);
        cpx_dvg = cpx_dv(gi);
        cpy_dvg = cpy_dv(gi);
        dx_dvg = dx_d(gi);
        
        gi = find(abs(zscore_mod(cpax_dv)) <= 3.5 & ...
                  abs(zscore_mod(cpay_dv)) <= 3.5);
        cpax_dvg = cpax_dv(gi);
        cpay_dvg = cpay_dv(gi);
        dxa_dvg = dx_d(gi);
        
        % mean and standard deviation of phase speed estimates        
        cpx(ii) = mean(cpx_dvg);
        cpy(ii) = mean(cpy_dvg);
        cpx_std(ii) = std(cpx_dvg);
        cpy_std(ii) = std(cpy_dvg);        
        
        cpax(ii) = mean(cpax_dvg);
        cpay(ii) = mean(cpay_dvg);
        cpax_std(ii) = std(cpax_dvg);
        cpay_std(ii) = std(cpay_dvg);   
        
        whm_event(ii) = nanmean(whm_tide_d);
        uhmc_event(ii) = nanmean(uhmc_tide_d);        
        
        cpa_event(ii) = cpa;
        cp_event(ii) = cp;    
    
        daten_events(ii) = event_daten(jj);
        
        t3_events(:,ii) = tcal3(di);
        t5_events(:,ii) = t5(di);
        tcrn_events(:,ii) = tcrn(di);
        tend_events(:,ii) = tend(di);
        tc_events(:,:,ii) = tc(di,:);
        tsgc_events(:,:,ii) = tsgc(di,:);
        tdts_events(:,:,ii) = tempC(:,di)'; 
        wc_events(:,:,ii) = C.M.evm(:,Ci)' + i*C.M.nvm(:,Ci)'; 
        wcdt_events(:,:,ii) = C.M.evm(:,Ci)' + i*C.M.nvm(:,Ci)' - C.utide(:,Ci)' - i*C.vtide(:,Ci)'; 
        
        %%% Make plots of three individual events
        zt_all = [zsC_isle 15];
        if plot_event
            set(gcf, 'PaperSize', [6.0 9.5]);
            set(gcf, 'PaperPosition', [0 0 6.5 5])   
            
            event_match = ii == pii;
            if any(event_match)
                col = find(event_match);
                if col == 1
                    labels = ['a)';'b)';'c)';'d)'];
                elseif col == 2
                    labels = ['e)';'f)';'g)';'h)'];
                else 
                    labels = ['i)';'j)';'k)';'l)'];
                end
                    
                subplot(7,3,col)
                ddi = di(1):di(end)+120;
                plot(datetime(ddi),t5(ddi),'r','linewidth',1.5);
                hold on
                plot(datetime(ddi),tcrn(ddi),'k','linewidth',1.5);
                plot(datetime(ddi),tend(ddi),'b','linewidth',1.5);
                ms = 4;
                plot(datetime(di(i5)),t5(di(i5)),'ro','markerfacecolor','r','markersize',ms);
                plot(datetime(di(ic)),tcrn(di(ic)),'ko','markerfacecolor','k','markersize',ms);
                plot(datetime(di(ie)),tend(di(ie)),'bo','markerfacecolor','b','markersize',ms);
                title({datestr(event_daten(jj))},'fontsize',10)
                hold off
                xlim([t1 t2]);
                xl = xlim;
                ymax = max([t5(ddi);tcrn(ddi);tend(ddi)]);
                ymin = min([t5(ddi);tcrn(ddi);tend(ddi)]);
                ylim([ymin-0.1 ymax+0.1])
                datetick('x','keeplimits')
                set(gca,'XTickLabel',[])
                if col == 1
                    ylabel('[deg C]')
                    fs = 8;
                    text(event_daten(jj)+0.08,19.75,'onshore','fontsize',fs,'color','b')
                    text(event_daten(jj)+0.08,19.55,'corner','fontsize',fs,'color','k')
                    text(event_daten(jj)+0.08,19.35,'east','fontsize',fs,'color','r')
                end
                set(gca,'tickdir','out')
                xl = double(xlim);
                yl = double(ylim);
                text(xl(1)-0.18*diff(xl),yl(2)+0.22*diff(yl),labels(1,:))
                if col == 3
                    text(xl(2)+0.03,yl(2),'DTS T')
                end

                ax2 = subplot(7,3,[3,6]+col);
                pcolor(datetime(di),distance/1000,tempC(:,di)-tcrn(di(ic))) 
                shading flat
                datetick('x')
                colormap(ax2,jet)
                xlim(xl)
                hold on
                plot(xl,[distance(zi5)/1000 distance(zi5)/1000],'k--')
                plot(xl,[distance(zicr)/1000 distance(zicr)/1000],'k--')
                plot(xl,[distance(zie)/1000 distance(zie)/1000],'k--')
                hold off
                caxis([-0.9,0.9])
                set(gca,'XTickLabel',[])
                if col == 1
                    ylabel({'distance','d [km]'})
                end
                if col == 3
                    cbar = colorbar('east');
                    pos = get(cbar,'position');
                    shift = 0.1;
                    shifty = -0.03;
                    sqz = 0.5;
                    sqzy = 0.75;
                    pos(1) = pos(1)+shift;
                    pos(2) = pos(2)+shifty;
                    pos(3) = pos(3)*sqz;
                    pos(4) = pos(4)*sqzy;
                    set(cbar,'position',pos)
                end
                set(gca,'tickdir','out')
                xl = double(xlim);
                yl = double(ylim);
                text(xl(1)-0.18*diff(xl),yl(2)+0.02*diff(yl),labels(2,:))
                if col == 3
                    text(xl(2)+0.03,yl(2),{'DTS T','anomaly','[deg C]'},'verticalalignment','top')
                end
                
                ax3 = subplot(7,3,[9,12]+col);
                pcolor(C.M.mtime,C.M.z+zoff,real(wc));
                shading flat
                caxiscen;
                cax = caxis;
                caxis([-.3 .3])
                hold on
                tc_all = [tc(di,:) tsgc(di)];
                contour(datetime(di),15-zt_all,tc_all',[0:.25:100],'color',[0.5 0.5 0.5])
                [c,h] = contour(datetime(di),15-zt_all,tc_all',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',1.5);
                clabel(c,h,'color',[0.3 0.3 0.3],'labelspacing',1000)
                hold off
                colormap(ax3,redblue)
                xlim(xl)
                ylim([0 15])
                datetick('x','keeplimits')
                set(gca,'layer','top')
                set(gca,'XTickLabel',[])
                if col == 1
                    ylabel('[mab]')
                end
                if col == 3
                    cbar = colorbar('east');
                    pos = get(cbar,'position');
                    pos(1) = pos(1)+shift;
                    pos(2) = pos(2)+shifty;
                    pos(3) = pos(3)*sqz;
                    pos(4) = pos(4)*sqzy;                    
                    set(cbar,'position',pos)
                end
                set(gca,'tickdir','out')
                xl = double(xlim);
                yl = double(ylim);
                text(xl(1)-0.28*diff(xl),yl(2)+0.02*diff(yl),labels(3,:))
                if col == 3
                    text(xl(2)+0.03,yl(2),{'east','velocity','[m/s]'},'verticalalignment','top')
                end
                
                ax4 = subplot(7,3,[15,18]+col);
                pcolor(C.M.mtime,C.M.z+zoff,imag(wc));
                shading flat
                caxiscen;
                cax = caxis;
                caxis([-.15 .15])
                hold on
                tc_all = [tc(di,:) tsgc(di)];
                zt_all = [zsC_isle 15];
                contour(datetime(di),15-zt_all,tc_all',[0:.25:100],'color',[0.5 0.5 0.5],'linewidth',0.5)
                [c,h] = contour(datetime(di),15-zt_all,tc_all',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',1.5);
                clabel(c,h,'color',[0.3 0.3 0.3],'labelspacing',1000)
                hold off
                colormap(ax4,redblue)
                xlim(xl)
                ylim([0 15])
                datetick('x','keeplimits')
                set(gca,'layer','top')
                if col == 1
                    ylabel('[mab]')
                end
                if col == 3
                    cbar = colorbar('east');
                    pos = get(cbar,'position');
                    pos(1) = pos(1)+shift;
                    pos(2) = pos(2)+shifty;
                    pos(3) = pos(3)*sqz;
                    pos(4) = pos(4)*sqzy;                    
                    set(cbar,'position',pos)
                end
                set(gca,'tickdir','out')
                xl = double(xlim);
                yl = double(ylim);
                text(xl(1)-0.28*diff(xl),yl(2)+0.02*diff(yl),labels(4,:))  
                if col == 3
                    text(xl(2)+0.03,yl(2),{'north','velocity','[m/s]'},'verticalalignment','top')
                end
            end

        end
        ii = ii+1;     
    end
end

print('-cmyk','-dpng',['../figures/fig_temp_velocity_c_events'])

%% Scatter plots of phase velocity
figure

subplot(2,2,2)
set(gcf, 'PaperSize', [7.0 2.7]);
set(gcf, 'PaperPosition', [0 -3.1 7.0 6.0])
h = errorbarxy(cpax,cpay,cpax_std,cpay_std,{'ko', 'k', 'k'});
axis equal
hz = zeroline('xy');
set(h.hMain,'MarkerFaceColor','r')
xlabel('c^\prime_{p}^x [m/s]')
ylabel('c^\prime_{p}^y [m/s]')
title('phase velocity relative to mean flow')
text(cpax(1),cpay(1),'A ',...
    'verticalalignment','top',...
    'horizontalalignment','right')
text(cpax(5),cpay(5),'B ',...
    'verticalalignment','top',...
    'horizontalalignment','right')
text(cpax(15),cpay(15)-0.003,'C ',...
    'verticalalignment','top',...
    'horizontalalignment','left')
xl = xlim;
yl = ylim;
text(xl(1)+0.02*diff(xl),yl(2)-0.05*diff(yl),'b)')

subplot(221)
h = errorbarxy(cpx,cpy,cpx_std,cpy_std,{'ko', 'k', 'k'});
axis equal
hz = zeroline('xy');
set(h.hMain,'MarkerFaceColor','r')
xlabel('c_{p}^x [m/s]')
ylabel('c_{p}^y [m/s]')
title('absolute phase velocity')
text(cpx(1),cpy(1),'A ',...
    'verticalalignment','top',...
    'horizontalalignment','right')
text(cpx(5),cpy(5),'B ',...
    'verticalalignment','top',...
    'horizontalalignment','right')
text(cpx(15),cpy(15)-0.003,'C ',...
    'verticalalignment','top',...
    'horizontalalignment','left')
xl = xlim;
yl = ylim;
text(xl(1)+0.02*diff(xl),yl(2)-0.05*diff(yl),'a)')

%%% Print phase velocity information for three events plotted
evi = [1,5,15];

for ii = 1:length(evi)
    disp(datestr(daten_events(evi(ii))))
    disp('-----------------------')
    disp(['cp_x = ' num2str(cpx(evi(ii)),3)])
    disp(['cp_y = ' num2str(cpy(evi(ii)),3)])
    disp(['cp_x (adj) = ' num2str(cpax(evi(ii)),3)])
    disp(['cp_y (adj) = ' num2str(cpay(evi(ii)),3)])
end

print('-dpdf','../figures/fig_cp_scatter')

%%
t0i = find(mitime_events==0);
dt0i = find(datetime_events==0);
zsi = 18; % surface adcp bin

% bottom velocity at front arrival
wb_events = squeeze(wc_events(t0i,1,:));

% 
g = 9.8;
H = 15; % water depth
alpha = 2.629e-4; % thermal expansion coeff (from T_dens_regression.m)
rhointcp = 1.0277e+03; % y-intercept in linear EOS
rho0 = 1025;

% stratification estimate
Tc_events = [tc_events tsgc_events];
lagi = 6; % -6 is a 1 hour lead
zti = 3;
deltaTc = squeeze(Tc_events(dt0i-lagi,end,:)-Tc_events(dt0i-lagi,1,:));
deltazc = zt_all(end-1)-zt_all(end);
dTdzc = deltaTc./deltazc;
deltarhoc = deltaTc*rho0*alpha;
drhodzc = deltarhoc./deltazc;
N = sqrt((g/1025)*drhodzc);

%temporal temperature difference
T1 = tsgc_events(dt0i-6,:); % before
T2 = tsgc_events(dt0i+6,:); % after
Tdiff = T2-T1;
rhodiff = -alpha*rho0*Tdiff;

% potential energy calculations 
ht_all = [15,15-zt_all];
dz = diff(ht_all);
ht_all_2d = repmat(ht_all',[1,size(Tc_events,3)]);
rho_t1_tmp = squeeze(rhointcp-alpha*rho0*Tc_events(dt0i-6,:,:));
rho_t2_tmp = squeeze(rhointcp-alpha*rho0*Tc_events(dt0i+6,:,:));
rho_t1 = [rho_t1_tmp(1,:); rho_t1_tmp];
rho_t2 = [rho_t2_tmp(1,:); rho_t2_tmp];
PE_t1 = trapz(ht_all,rho_t1*g.*ht_all_2d);
PE_t2 = trapz(ht_all,rho_t2*g.*ht_all_2d);
delta_PE = PE_t2-PE_t1;

%horizontal temperature difference
tdts2 = tsgc_events(dt0i,:);
tdts1 = t5_events(dt0i,:);

tdts1 = squeeze(tdts_events(dt0i,zic-20,:))';
tdts2 = squeeze(tdts_events(dt0i,zic+20,:))';

Tdiffh = tdts2-tdts1;
rhodiffh = alpha*rho0*Tdiffh;

deltaT2 = squeeze(Tc_events(dt0i+lagi,2:end,:)-Tc_events(dt0i+lagi,1:end-1,:));
deltaz = zt_all(2:end)-zt_all(1:end-1);
deltaz2 = repmat(deltaz,[size(deltaT2,2),1]);
dTdz2 = -deltaT2'./deltaz2;
[~,si] = max(dTdz2,[],2);

h = 7.5; % h = 1/2 d, energy conserving case (Benjamin 1968)

c0 = sqrt(2*g*h.*rhodiff/rho0);
S = -rhodiffh./deltarhoc';

cpmag = sqrt(cpx.^2+cpy.^2);
cpamag = sqrt(cpax.^2+cpay.^2);

phi = angle(cpx+i*cpy);
um = real(whm_event.*exp(-i*phi));
cpa2 = cpmag - um;
cpamag2 = abs(cpa2);
cpax2 = cpamag2.*cos(phi);
cpay2 = cpamag2.*sin(phi);

cwave = N*15/3.14;

% linear regression analysis
gi = isfinite(cpamag+rhodiff);
[r,p] = corrcoef(cpmag(gi),c0(gi));
m = polyfit(c0(gi),cpmag(gi),1);

gi = isfinite(cpamag+rhodiff);
[ra,pa] = corrcoef(cpmag(gi),c0(gi));
ma = polyfit(c0(gi),cpamag(gi),1);

gi = isfinite(cpamag+rhodiff);
[ry,py] = corrcoef(cpy(gi),c0(gi));
my = polyfit(c0(gi),cpy(gi),1);

gi = isfinite(cpamag+rhodiff);
[ra,pa] = corrcoef(cpamag(gi),c0(gi));
ma = polyfit(c0(gi),cpamag(gi),1);

gi = isfinite(cpamag+rhodiff);
[ray,pay] = corrcoef(cpay(gi),c0(gi));
may = polyfit(c0(gi),cpay(gi),1);


figure()
set(gcf,'papersize',[6,3])
set(gcf,'paperposition',[-0.2 -2.6 6.5 5.5])


lw = 2; % line width
lc = [0.5 0.5 0.5]; % line color
tc = lc; % text color

subplot(221)
scatter(c0,cpy,30,cwave,'filled','k')
hold on 
plot(c0,c0*my(1)+my(2),'color',lc,'linewidth',lw)
box on
xlabel('gravity current c_{gc} [m/s]')
ylabel('observed c_p^y [m/s]')
title({'phase speed','onshore component'})
xl = xlim;
yl = ylim;
text(0.35*diff(xl),yl(2)+0.06*diff(yl),'a)')
text(0.95*diff(xl),yl(1)+0.33*diff(yl),['m = ', num2str(my(1),2)],'color',tc)
text(0.95*diff(xl),yl(1)+0.25*diff(yl),['b = ', num2str(my(2),1)],'color',tc)
text(0.95*diff(xl),yl(1)+0.14*diff(yl),['r = ', num2str(ry(2),2)],'color',tc)
text(0.95*diff(xl),yl(1)+0.06*diff(yl),['p = ', num2str(py(2),1)],'color',tc)

subplot(222)
scatter(c0,cpay,30,cwave,'filled','k')
hold on 
plot(c0,c0*may(1)+may(2),'color',lc,'linewidth',lw)
xlabel('gravity current c_{gc} [m/s]')
ylabel('observed c''_p^{ y} [m/s]')
title({'adjusted phase speed','onshore component'})
box on
xl = xlim;
yl = ylim;
text(0.36*diff(xl),yl(2)-0.07*diff(yl),'b)')
text(0.95*diff(xl),yl(1)+0.23*diff(yl),['m = ', num2str(may(1),2)],'color',tc)
text(0.95*diff(xl),yl(1)+0.16*diff(yl),['b = ', num2str(may(2),1)],'color',tc)
text(0.95*diff(xl),yl(1)+0.05*diff(yl),['r = ', num2str(ray(2),2)],'color',tc)
text(0.95*diff(xl),yl(1)-0.02*diff(yl),['p = ', num2str(pay(2),1)],'color',tc)

print('-dpdf','../figures/fig_cp_compare')

%%
figure
subplot(221)
scatter(c0,cpamag2,30,cwave,'filled','k')
hold on
plot(c0,c0*m(1)+m(2),'color',lc,'linewidth',lw)
box on
xlabel('gravity current c_{gc} [m/s]')
ylabel('observed c [m/s]')
title({'phase speed','magnitude'})
xl = xlim;
yl = ylim;
text(0.35*diff(xl),yl(2)-0.07*diff(yl),'a)')
text(0.95*diff(xl),0.13*diff(yl),['m = ', num2str(m(1),3)],'color',tc)
text(0.95*diff(xl),0.06*diff(yl),['b = ', num2str(m(2),1)],'color',tc)
text(0.95*diff(xl),-0.05*diff(yl),['r = ', num2str(r(2),2)],'color',tc)
text(0.95*diff(xl),-0.12*diff(yl),['p = ', num2str(p(2),1)],'color',tc)
