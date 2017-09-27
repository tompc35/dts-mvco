clear all

plot_event = 1;

load_dts_isle_data
detide_c;

nfilt = 3;

% interpolate seagauge temp at C to DTS time base
tsgc = interp1(mdaysg,wtsg(:,3),datetime);
[eventi,event_daten] = get_event_indices_dTdt(boxfilt(tsgc,nfilt),datetime,-0.5);

%%

% 5 - calibration point 5
% c - corner of DTS cable
% e - end of DTS cable


ti = 1:length(eventi);

th = interp1(mday_isle,wtH_isle,datetime);
tc = interp1(mday_isle,wtC_isle,datetime);
te = interp1(mday_isle,wtE_isle,datetime);
tI = interp1(mday_isle,wtI_isle,datetime);
tG = interp1(mday_isle,wtG_isle,datetime);
tnode = interp1(mday_node,wtb_node,datetime);

close all
clear c_*
clear cpx cpy cpax cpay
clear *_events

ii = 1;
for kk = 1:length(ti)
%for kk = 1:1
    
    jj = ti(kk);
    t1 = event_daten(jj)-0.25;
    t2 = event_daten(jj)+0.25;
    di = find(datetime>=t1 & datetime<=t2);
    hi = find(H.ttime>=t1 & H.ttime<t2);
    Ii = find(I.M.mtime>=t1 & I.M.mtime<t2);
    Ci = find(C.M.mtime>=t1 & C.M.mtime<t2);
    tim = find(mday_isle>=t1 & mday_isle<=t2);
    
    if ii == 1
        datetime_events = datetime(di)-event_daten(jj);
        mday_isle_events = mday_isle(tim)-event_daten(jj);
        ttime_events = H.ttime(hi)-event_daten(jj);
        mitime_events = I.M.mtime(Ii)-event_daten(jj);
    end
    
    % use only events where data is good (check corner)
    if ~isnan(tempC(zic,di))
    
        for zii = 1:40
            zjj = zii;
            %zii = 10;
            % modify points
            zicr = zic; %  corner
            zi5 = zic-zii;
            zie = zic+zjj;

            latc = lat_dts(zicr);
            lonc = lon_dts(zicr);

            late = lat_dts(zie);
            lone = lon_dts(zie);

            % geometry of mooring array
            xc = 0; yc = 0;
            [x5,y5] = latlon2xy(lat_dts(zi5),lon_dts(zi5),latc,lonc);
            [xe,ye] = latlon2xy(late,lone,latc,lonc);

            % angle of deviation of c-e line from x (eastward) axis
            theta = atan2(ye,xe)*180/pi;

            % length of triangle segments
            length_5c = sqrt((x5-xc)^2 + (y5-yc)^2);
            length_ce = sqrt((xc-xe)^2 + (yc-ye)^2);
            length_e5 = sqrt((xe-x5)^2 + (ye-y5)^2);

            % angles of triangle (law of cosines)
            ang_5c = acos((-length_5c^2 + length_ce^2 + length_e5^2) ...
                            /(2*length_ce*length_e5))*180/pi;
            ang_ce = acos((-length_ce^2 + length_e5^2 + length_5c^2) ...
                            /(2*length_e5*length_5c))*180/pi;      
            ang_e5 = acos((-length_e5^2 + length_5c^2 + length_ce^2) ...
                            /(2*length_5c*length_ce))*180/pi;  

            tcrn = tempC(zicr,:);
            tend = tempC(zie,:);
            %t5 = tcal5';t5 
            t5 = tempC(zi5,:);

            tcrn = boxfilt(tcrn,nfilt);
            tend = boxfilt(tend,nfilt);
            t5 = boxfilt(t5,nfilt);

            dt = (datetime(2)-datetime(1))*86400;

            dTdtcrn = nan(size(tcrn));
            dTdtcrn(2:end-1) = 0.5*(tcrn(3:end)-tcrn(1:end-2))/dt;

            dTdtend = nan(size(tend));
            dTdtend(2:end-1) = 0.5*(tend(3:end)-tend(1:end-2))/dt;

            dTdt5 = nan(size(t5));
            dTdt5(2:end-1) = 0.5*(t5(3:end)-t5(1:end-2))/dt;


            [~,i5] = min(dTdt5(di));
            [~,ic] = min(dTdtcrn(di));
            [~,ie] = min(dTdtend(di));

            tstart = datetime(di(min([i5 ic ie])));
            tf = datetime(di(max([i5 ic ie])));

            cp_5c = 1000*length_5c/((i5-ic)*dt);
            cp_ce = -1000*length_ce/((ic-ie)*dt);
            cp_e5 = 1000*length_e5/((ie-i5)*dt);    

            phi = atan((cp_ce/cp_5c)*cscd(ang_e5) - cotd(ang_e5));
            phi2 = atan(cotd(ang_5c) - (cp_ce/cp_e5)*cscd(ang_5c));

            %disp(datestr(datetime(di(1))))
            %disp(['phi_1 = ' num2str(phi*180/pi)])
            %disp(['phi_2 = ' num2str(phi2*180/pi)])

            % magnitude of phase velocity
            cp1 = cp_ce*cos(phi);
            cp2 = cp_5c*cos(ang_e5*pi/180 - phi);
            cp3 = cp_e5*cos(phi+ang_5c*pi/180);

            if isfinite(cp1)
                cp = cp1;
            elseif isfinite(cp2)
                cp = cp2;
            elseif isfinite(cp3)
                cp = cp3;
            else 
                cp = NaN;
            end

            if cp < 0
               cp = -cp;
               phi = phi - pi;
            end

            % angle of propagation relative to x-axis
            phixy = theta*pi/180-phi; 

            % mean tidal velocity in direction of phase propagation
             cmi = find(C.M.mtime>=tstart & C.M.mtime<tf);

            zoff = 0.75; % difference between bottom and ADCP height

            %uhm_tide = mean(depthavg(C.utide(1:17,cmi),15-(C.M.z(1:17)+zoff),15),2);
            %vhm_tide = mean(depthavg(C.vtide(1:17,cmi),15-(C.M.z(1:17)+zoff),15),2);

            uhm_tide = mean(depthavg(C.M.evm(1:17,cmi),15-(C.M.z(1:17)+zoff),15),2);
            vhm_tide = mean(depthavg(C.M.nvm(1:17,cmi),15-(C.M.z(1:17)+zoff),15),2);

            whm_tide = uhm_tide+i*vhm_tide;

            whmc_tide = whm_tide*exp(-i*phixy);
            uhmc_tide = real(whmc_tide);   

            wc = C.M.evm + i*C.M.nvm;
            wcc = wc*exp(-i*phixy);


            wcdt = C.M.evm + i*C.M.nvm - C.utide - i*C.vtide;
            wccdt = wcdt*exp(-i*phixy);

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
        end
        
        %n = 20;
        
        m = 1;
        n = 40;
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
        
        %ch_crit = norminv(1-1/(4*length(cpax_dv)),0,1);
        
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
        
        if ii==1
            
            figure(1)
            set(gcf, 'PaperSize', [7.0 7]);
            set(gcf, 'PaperPosition', [0 0 7 7])   

            subplot(221)
            h1 = plot(dxa_dvg,cpax_dvg,'-','color',[.7 .7 .7],'linewidth',2);
            hold on
            h2 = plot(dxa_dvg,cpay_dvg,'k-','linewidth',2);
            xlabel('\Deltad in onshore direction [m]')
            ylabel(['[m/s]'])
            leg = legend([h1(1);h2(1)],'c^\prime_{p}^x','c^\prime_{p}^y','location','northwest');

%             subplot(222)
%             plot(dy_d,cpay_d,'k.')
%             hold on
%             plot(dy_d,cpax_d,'.','color',[.7 .7 .7])
%             xlabel('\Deltad in east direction [m]')
%             ylabel('[m/s]')
% 
%             subplot(2,2,[3:4])
%             h = histogram(cpay_dv,40,'normalization','count');
%             hold on
%             hh = histogram(cpax_dv,40,'normalization','count');
%             h.FaceColor = [0 0 0];
%             hh.FaceColor = [1 1 1];
%             title('histogram for \Deltad \leq 75m')
%             legh = legend([hh(1);h(1)],'c^\prime_{p}^x','c^\prime_{p}^y','location','northwest');
%             ylabel('count')
%             xlabel('[m/s]')

            title(datestr(event_daten(jj)))
            hold off


            %print('-dpdf',['figures_paper/fig_event_cp_stats_' datestr(datetime(di(1)),'YYYYMMDD')])

        end
        
        
        %%
        
%         cpx(ii) = cp*cos(phixy);
%         cpy(ii) = cp*sin(phixy);
%         
%         cpax(ii) = cpa*cos(phixy);
%         cpay(ii) = cpa*sin(phixy);
        
        cpa_event(ii) = cpa;
        cp_event(ii) = cp;    
        %phixy_event(ii) = phixy;
    
        daten_events(ii) = event_daten(jj);
        t3_events(:,ii) = tcal3(di);
        t5_events(:,ii) = t5(di);
        tcrn_events(:,ii) = tcrn(di);
        tend_events(:,ii) = tend(di);
        tnode_events(:,ii) = tnode(di);
        th_events(:,:,ii) = th(di,:);
        tc_events(:,:,ii) = tc(di,:);
        tsgc_events(:,:,ii) = tsgc(di,:);
        te_events(:,:,ii) = te(di,:);
        tI_events(:,:,ii) = tI(di,:);
        tG_events(:,:,ii) = tG(di,:);
        tdts_events(:,:,ii) = tempC(:,di)';
        wi_events(:,:,ii) = I.M.evm(:,Ii)' + i*I.M.nvm(:,Ii)';   
        wg_events(:,:,ii) = G.M.evm(:,Ii)' + i*G.M.nvm(:,Ii)'; 
        wc_events(:,:,ii) = C.M.evm(:,Ii)' + i*C.M.nvm(:,Ii)'; 
        wcdt_events(:,:,ii) = C.M.evm(:,Ii)' + i*C.M.nvm(:,Ii)' - C.utide(:,Ii)' - i*C.vtide(:,Ii)'; 
        intensc_events(:,:,ii) = 0.25*(C.M.intens1(:,Ii)' + C.M.intens2(:,Ii)'+C.M.intens3(:,Ii)' + C.M.intens4(:,Ii)');   
        
        %%
        zt_all = [zsC_isle 15];
        if plot_event
            figure
            set(gcf, 'PaperSize', [8.5 6.0]);
            set(gcf, 'PaperPosition', [0 0 5 7.0])

            subplot(711)
            ddi = di(1):di(end)+120;
            plot(datetime(ddi),t5(ddi),'r','linewidth',2);
            hold on
            plot(datetime(ddi),tcrn(ddi),'k','linewidth',2);
            plot(datetime(ddi),tend(ddi),'b','linewidth',2);
            plot(datetime(di(i5)),t5(di(i5)),'ro','markerfacecolor','r');
            plot(datetime(di(ic)),tcrn(di(ic)),'ko','markerfacecolor','k');
            plot(datetime(di(ie)),tend(di(ie)),'bo','markerfacecolor','b');
            title({datestr(event_daten(jj)),...
                        'DTS temperature at three points'})
            leg = legend('east','corner','onshore','location','northeast');
            set(leg,'fontsize',8)
            hold off
            xlim([t1 t2+4/24]);
            xl = xlim;
            ymax = max([t5(ddi);tcrn(ddi);tend(ddi)]);
            ymin = min([t5(ddi);tcrn(ddi);tend(ddi)]);
            ylim([ymin-0.1 ymax+0.1])
            datetick('x','keeplimits')
            ylabel('[deg C]')
            set(gca,'XTickLabel',[])


            ax2 = subplot(7,1,[2:3]);
            pcolor(datetime(di),distance/1000,tempC(:,di)) 
            shading flat
            datetick('x')
            colormap(ax2,jet)
            colorbar('east')
            xlim(xl)
            hold on
            plot(xl,[distance(zi5)/1000 distance(zi5)/1000],'k--')
            plot(xl,[distance(zicr)/1000 distance(zicr)/1000],'k--')
            plot(xl,[distance(zie)/1000 distance(zie)/1000],'k--')
            hold off
            ylabel('distance, d [km]')
            title('DTS temperature [^oC]')
            set(gca,'XTickLabel',[])

            ax3 = subplot(7,1,[4:5]);
            pcolor(C.M.mtime,C.M.z+zoff,real(wc));
            %pcolor(C.M.mtime,C.M.z,C.M.evm);
            shading flat
            caxiscen;
            cax = caxis;
            caxis([-.4 .4])
            hold on
            tc_all = [tc(di,:) tsgc(di)];
            contour(datetime(di),15-zt_all,tc_all',[0:.1:100],'color',[0.5 0.5 0.5])
            [c,h] = contour(datetime(di),15-zt_all,tc_all',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2);
            clabel(c,h,'color',[0.3 0.3 0.3],'labelspacing',1000)
            hold off
            colormap(ax3,redblue)
            colorbar('east')
            xlim(xl)
            ylabel('[mab]')
            ylim([0 15])
            title({'eastward velocity [m/s]'})
            datetick('x','keeplimits')
            set(gca,'layer','top')
            set(gca,'XTickLabel',[])

            ax4 = subplot(7,1,[6:7]);
            pcolor(C.M.mtime,C.M.z+zoff,imag(wc));
            %pcolor(C.M.mtime,C.M.z,C.M.nvm);
            shading flat
            caxiscen;
            cax = caxis;
            caxis([-.1 .1])
            hold on
            tc_all = [tc(di,:) tsgc(di)];
            zt_all = [zsC_isle 15];
            contour(datetime(di),15-zt_all,tc_all',[0:.1:100],'color',[0.5 0.5 0.5])
            [c,h] = contour(datetime(di),15-zt_all,tc_all',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2);
            clabel(c,h,'color',[0.3 0.3 0.3],'labelspacing',1000)
            hold off
            colormap(ax4,redblue)
            colorbar('east')
            xlim(xl)
            ylabel('[mab]')
            ylim([0 15])
            title({'northward velocity [m/s]'}) 
            datetick('x','keeplimits')
            set(gca,'layer','top')

            print('-cmyk','-dpng',['figures_paper/events/fig_temp_velocity_' datestr(datetime(di(1)),'yyyymmdd')])
            %print('-dpng',['figures/h_events_agu/h_event_propagation_profile_' datestr(datetime(di(1)),'mmdd_HHMM')])
        
%             figure
%             quiver(cpax(ii),cpay(ii))
%             hold on
%             quiver(cpx(ii),cpy(ii))
        end
        %%
        ii = ii+1;
        
    end
end

%%
figure
set(gcf, 'PaperSize', [8.5 6.0]);
set(gcf, 'PaperPosition', [0 0 5 7.0])
ax1 = subplot(311);
pcolor(datetime_events,distance/1000,squeeze(nanmean(tdts_events,3))')
yl = ylim;
xl = [-0.25 0.25];
xlim(xl);
shading flat
colormap(ax1,jet)
hold on
plot(xl,[distance(zicr) distance(zicr)]/1000,'k--')
colorbar('eastoutside')
ylabel('distance, d [km]')
box on
title('DTS temperature [^oC]')

subplot(312)
pcolor(mitime_events,C.M.z+0.75,real(squeeze(mean(wcdt_events,3)))'), shading flat
caxiscen;
cax = caxis;
hold on
contour(datetime_events, 15-zt_all, squeeze(mean([tc_events tsgc_events],3))',[0:.1:100],'color',[0.5 0.5 0.5])
[c,h] = contour(datetime_events, 15-zt_all, squeeze(mean([tc_events tsgc_events],3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
clabel(c,h)
ylim([0 15])
colorbar('eastoutside')
xlim(xl);
ylabel('[mab]')
title('de-tided eastward velocity [m/s]')


subplot(313)
pcolor(mitime_events,C.M.z+0.75,imag(squeeze(mean(wcdt_events,3)))'), shading flat
caxiscen;
cax = caxis;
hold on
contour(datetime_events, 15-zt_all, squeeze(mean([tc_events tsgc_events],3))',[0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events, 15-zt_all, squeeze(mean([tc_events tsgc_events],3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
ylim([0 15])
colorbar('eastoutside')
colormap(redblue)
xlim(xl);
ylabel('[mab]')
xlabel('time relative to event start [days]')
title('de-tided northward velocity [m/s]')

print('-cmyk','-dpng',['figures_paper/fig_temp_velocity_composite'])

%%
figure(1)
subplot(2,2,3)
set(gcf, 'PaperSize', [7.0 6.0]);
set(gcf, 'PaperPosition', [0 0 7.0 6.0])
h = errorbarxy(cpax,cpay,cpax_std,cpay_std,{'ko', 'k', 'k'});
axis equal
hz = zeroline('xy');
set(h.hMain,'MarkerFaceColor','r')
xlabel('c^\prime_{p}^x [m/s]')
ylabel('c^\prime_{p}^y [m/s]')
title('phase velocity relative to mean flow')
%print('-dpdf','figures_paper/fig_cp_scatter')


subplot(221)
h = errorbarxy(cpx,cpy,cpx_std,cpy_std,{'ko', 'k', 'k'});
axis equal
hz = zeroline('xy');
set(h.hMain,'MarkerFaceColor','r')
xlabel('c_{p}^x [m/s]')
ylabel('c_{p}^y [m/s]')
title('absolute phase velocity')
print('-dpdf','figures_paper/fig_cp_scatter')

%%
figure
subplot(211)
pcolor(mitime_events,C.M.z+0.75,real(squeeze(mean(wcdt_events,3)))'), shading flat
caxiscen;
cax = caxis;
hold on
contour(datetime_events, 15-zt_all, squeeze(mean([tc_events tsgc_events],3))',[0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events, 15-zt_all, squeeze(mean([tc_events tsgc_events],3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
ylim([0 15])
colorbar
colormap(redblue)

subplot(212)
pcolor(mitime_events,C.M.z+0.75,imag(squeeze(mean(wcdt_events,3)))'), shading flat
caxiscen;
cax = caxis;
hold on
contour(datetime_events, 15-zt_all, squeeze(mean([tc_events tsgc_events],3))',[0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events, 15-zt_all, squeeze(mean([tc_events tsgc_events],3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
ylim([0 15])
colorbar
colormap(redblue)

%Standard error
wbotse = std(real(squeeze(wcdt_events(:,1,:)))')./sum(isfinite(imag(squeeze(wcdt_events(:,1,:)))'));

tsgc_anom = squeeze(tsgc_events(:,1,:)) - repmat(nanmean(squeeze(tsgc_events(:,1,:))),[length(datetime_events),1])

tsgc_anom_se = std(tsgc_anom)./sum(isfinite(tsgc_anom))

%%
figure
subplot(211)
plot(mitime_events,real(squeeze(wcdt_events(:,1,:)))'), shading flat

subplot(212)
plot(mitime_events,imag(squeeze(wcdt_events(:,1,:)))'), shading flat

% figure
% compass(cpx,cpy,'r')
% % hold on
% % compass(cpax,cpay,'r')
% 
% figure
% compass(cpax,cpay,'b')


%%
t0i = find(mitime_events==0);
dt0i = find(datetime_events==0);
zsi = 18; % surface adcp bin

% bottom velocity at front arrival
wb_events = squeeze(wcdt_events(t0i,1,:));

% 
g = 9.8;
H = 15; % water depth
alpha = 2.629e-4; % thermal expansion coeff (from T_dens_regression.m)
rho0 = 1025;

% stratification estimate
Tc_events = [tc_events tsgc_events];
lagi = -24; % -12 is a 1 hour lead
%zti = length(zt_all);
zti = 3;
deltaTc = squeeze(Tc_events(dt0i-lagi,end,:)-Tc_events(dt0i-lagi,1,:));
deltazc = zt_all(1)-zt_all(end);
dTdzc = deltaTc./deltazc;
deltarhoc = deltaTc*rho0*alpha;
drhodzc = deltaTc*rho0*alpha./deltazc;

%temporal temperature difference
T1 = tsgc_events(dt0i+6,:);
T2 = tsgc_events(dt0i-6,:);
Tdiff = T2-T1;
rhodiff = alpha*rho0*Tdiff;

%horizontal temperature difference
tdts2 = tsgc_events(dt0i,:);
%tdts1 = squeeze(tdts_events(dt0i,zie-10,:))';
tdts1 = t5_events(dt0i,:);
Tdiffh = tdts2-tdts1;
rhodiffh = alpha*rho0*Tdiffh;

c0 = sqrt(abs(g*rhodiffh/rho0));
S = -rhodiff./deltarhoc';!

cpmag = sqrt(cpx.^2+cpy.^2);
cpamag = sqrt(cpax.^2+cpay.^2);

figure
scatter(cpmag,imag(wb_events),30,deltarhoc,'filled');

figure
plot(cpmag,drhodzc,'.');

figure
scatter(cpmag,c0,30,log(S),'filled');
gi = isfinite(cpmag+rhodiffh);
[r,p] = corrcoef(cpmag(gi),sqrt(c0(gi)));

figure
plot(cpamag,sqrt(drhodzc*9.8/1025)*15,'.');