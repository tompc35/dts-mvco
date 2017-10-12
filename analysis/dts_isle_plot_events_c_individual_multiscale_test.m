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

% indices of events to plot (three)
pii = [1,5,15];

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
    
        for zi = 1:40
            %zii = zi+40;
            zii = zi;
            zjj = zii;
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

            dx_d(zi) = (distance(end)-distance(end-1))*(zii);
            dy_d(zi) = (distance(end)-distance(end-1))*(zjj);

            ur_d(zi) = mean(uhmc_tide);
            u_d(zi) = mean(real(whm_tide));
            v_d(zi) = mean(imag(whm_tide));

            cp_d(zi) = cp;
            cpa_d(zi) = cpa;
            phixy_d(zi) = phixy;

            cpx_d(zi) = cp*cos(phixy);
            cpy_d(zi) = cp*sin(phixy);

            cpax_d(zi) = cpa*cos(phixy);
            cpay_d(zi) = cpa*sin(phixy);
        end
        
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

        %%
        
        cpa_event(ii) = cpa;
        cp_event(ii) = cp;    
    
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
                    squeeze = 0.5;
                    squeezey = 0.75;
                    pos(1) = pos(1)+shift;
                    pos(2) = pos(2)+shifty;
                    pos(3) = pos(3)*squeeze;
                    pos(4) = pos(4)*squeezey;
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
                contour(datetime(di),15-zt_all,tc_all',[0:.1:100],'color',[0.5 0.5 0.5])
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
                    pos(3) = pos(3)*squeeze;
                    pos(4) = pos(4)*squeezey;                    
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
                caxis([-.3 .3])
                hold on
                tc_all = [tc(di,:) tsgc(di)];
                zt_all = [zsC_isle 15];
                contour(datetime(di),15-zt_all,tc_all',[0:.1:100],'color',[0.5 0.5 0.5],'linewidth',0.5)
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
                    pos(3) = pos(3)*squeeze;
                    pos(4) = pos(4)*squeezey;                    
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

print('-cmyk','-dpng',['figures_paper/fig_temp_velocity_c_events'])



%%
figure
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

subplot(221)
h = errorbarxy(cpx,cpy,cpx_std,cpy_std,{'ko', 'k', 'k'});
axis equal
hz = zeroline('xy');
set(h.hMain,'MarkerFaceColor','r')
xlabel('c_{p}^x [m/s]')
ylabel('c_{p}^y [m/s]')
title('absolute phase velocity')
print('-dpdf','figures_paper/fig_cp_scatter')