clear all

load_dts_isle_data
%[eventi,event_daten] = get_event_indices(tcal3,tcal4,datetime,2);
nfilt = 3;
[eventi,event_daten] = get_event_indices_dTdt(boxfilt(tcal3,nfilt),datetime,-0.5);
%[eventi,event_daten] = get_event_indices_dTdt(boxfilt(tcal3,nfilt),datetime,-1.0);
detide_g;
detide_i;

%%

% find where H time series starts
% tmin=min(H.ttime(find(isfinite(sum(H.uu)))));
% ti = find(event_daten > tmin);
ti = 1:length(eventi);

th = interp1(mday_isle,wtH_isle,datetime);
tc = interp1(mday_isle,wtC_isle,datetime);
tG = interp1(mday_isle,wtG_isle,datetime);
tnode = interp1(mday_node,wtb_node,datetime);

% combine temperature mooring and quadpod at site E
tem = interp1(mday_isle,wtE_isle,datetime);
teq = fliplr(interp1(Equad.mtime,Equad.aT,datetime));
te = [tem, teq];
zsE_isle = [zsE_isle, 15-fliplr(Equad.mab)];

% combine temperature mooring and ADCP temp at site I
tIm = interp1(mday_isle,wtI_isle,datetime);
tIv = interp1(I.M.mtime,I.M.temp',datetime);
tI = [tIm,tIv];
zsI_isle = [zsI_isle, 20.63];

close all
for ii = 1:length(ti)
%for ii = 1:1
    jj = ti(ii);
    t1 = event_daten(jj)-0.3;
    t2 = event_daten(jj)+0.3;
    di = find(datetime>=t1 & datetime<=t2);
    hi = find(H.ttime>=t1 & H.ttime<t2);
    Ii = find(I.M.mtime>=t1 & I.M.mtime<t2);
    tim = find(mday_isle>=t1 & mday_isle<=t2);
    
    if ii == 1
        t3_events = nan([length(di) length(ti)]);
        th_events = nan([length(di) size(th,2) length(ti)]);
        tc_events = nan([length(di) size(tc,2) length(ti)]);
        te_events = nan([length(di) size(te,2) length(ti)]);
        tI_events = nan([length(di) size(tI,2) length(ti)]);
        tG_events = nan([length(di) size(tG,2) length(ti)]);
        tnode_events = nan([length(di) length(ti)]);
        tdts_events = nan([length(di) size(tempC,1) length(ti)]);
        wh_events = nan([length(hi) length(H.zz) length(ti)])+i*nan([length(hi) length(H.zz) length(ti)]);
        wh_detide_events = nan([length(hi) length(H.zz) length(ti)])+i*nan([length(hi) length(H.zz) length(ti)]);
        wi_events = nan([length(Ii) length(I.M.z) length(ti)]);
        wi_detide_events = nan([length(Ii) length(I.M.z) length(ti)]);   
        wg_events = nan([length(Ii) length(G.M.z) length(ti)]);
        wg_detide_events = nan([length(Ii) length(G.M.z) length(ti)]);       
        datetime_events = datetime(di)-event_daten(jj);
        mday_isle_events = mday_isle(tim)-event_daten(jj);
        ttime_events = H.ttime(hi)-event_daten(jj);
        mitime_events = I.M.mtime(Ii)-event_daten(jj);
    end
    
    t3_events(:,ii) = tcal3(di);
    tnode_events(:,ii) = tnode(di);
    th_events(:,:,ii) = th(di,:);
    tc_events(:,:,ii) = tc(di,:);
    te_events(:,:,ii) = te(di,:);
    tI_events(:,:,ii) = tI(di,:);
    tG_events(:,:,ii) = tG(di,:);
    tdts_events(:,:,ii) = tempC(:,di)';
    wh_events(:,:,ii) = H.uu(:,hi)' + i*H.vv(:,hi)';
    wh_detide_events(:,:,ii) = (H.uu(:,hi)-H.utide(:,hi))' +...
                                i*(H.vv(:,hi)-H.vtide(:,hi))';
    wi_events(:,:,ii) = I.M.evm(:,Ii)' + i*I.M.nvm(:,Ii)'; 
    wi_detide_events(:,:,ii) = I.M.evm(:,Ii)' + i*I.M.nvm(:,Ii)' - ...
                                (I.utide(:,Ii)+i*I.vtide(:,Ii))';
    wg_events(:,:,ii) = G.M.evm(:,Ii)' + i*G.M.nvm(:,Ii)'; 
    wg_detide_events(:,:,ii) = G.M.evm(:,Ii)' + i*G.M.nvm(:,Ii)' - ...
                                (G.utide(:,Ii)+i*G.vtide(:,Ii))'; 
end

t3_events_anom = t3_events-repmat(t3_events(round(length(di)/2),:),[length(di) 1]);
tnode_events_anom = tnode_events-repmat(tnode_events(round(length(di)/2),:),[length(di) 1]);


%%
% find events where velocity data exists at H
hidx = find(isfinite(squeeze(sum(wh_detide_events(15,:,:),2))));

% combine temperatures at H
t3_events_tmp(1,:,:) = t3_events;
t_events = [th_events permute(t3_events_tmp,[2 1 3])];
zt_events = [zsH_isle 15];



%%

figure
set(gcf, 'PaperSize', [7.0 7.0]);
set(gcf, 'PaperPosition', [0 0 6.8 6.2])

ax1 = subplot(15,2,[1:2:5])
pcolor(datetime_events*24,distance/1000,squeeze(nanmean(tdts_events,3))')
shading flat
colormap(ax1,jet)
colorbar
title('DTS temperature [^oC]')
ylabel('distance, d [km]')

subplot(15,2,[2:2:6])
contour(datetime_events*24,15-zsE_isle,squeeze(nanmean(te_events,3))',[0:.1:100],'color',[0.5 0.5 0.5])
hold on
contour(datetime_events*24,15-zsE_isle,squeeze(nanmean(te_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
title('E temperature [^oC]')
cbar = colorbar;
set(cbar,'visible','off')
ylabel('[mab]')
ylim([0,15])


subplot(15,2,[9:2:13])
pcolorjw(ttime_events*24,H.zz,squeeze(nanmean(real(wh_events(:,:,hidx)),3))')
%pcolorjw(ttime_events,H.zz,squeeze(nanmean(real(wh_events(:,:,hidx)),3))')
shading flat
colormap(redblue)
caxiscen;
hold on
contour(datetime_events*24,15-zt_events,squeeze(nanmean(t_events(:,:,hidx),3))', ...
            [0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events*24,15-zt_events,squeeze(nanmean(t_events(:,:,hidx),3))', ...
            [0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
colorbar
title('H - eastward velocity')
ylim([0,15])
caxis([-0.15,0.15])
ylabel('[mab]')

subplot(15,2,[10:2:14])
pcolorjw(ttime_events*24,H.zz,squeeze(nanmean(imag(wh_events(:,:,hidx)),3))')
%pcolorjw(ttime_events,H.zz,squeeze(nanmean(imag(wh_events(:,:,hidx)),3))')
shading flat
colormap(redblue)
caxiscen;
colorbar
hold on
contour(datetime_events*24,15-zt_events,squeeze(nanmean(t_events(:,:,hidx),3))', ...
            [0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events*24,15-zt_events,squeeze(nanmean(t_events(:,:,hidx),3))', ...
            [0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
title('H - northward velocity')
ylim([0,15])
caxis([-0.05,0.05])
ylabel('[mab]')

%%% G velocity %%%

zi = 1:15;
subplot(15,2,[17:2:19])
pcolorjw(mitime_events*24,G.M.z(zi),squeeze(nanmean(real(wg_events(:,zi,:)),3))')
%pcolorjw(mitime_events,G.M.z(zi),squeeze(nanmean(real(wg_events(:,zi,:)),3))')
shading flat
colormap(redblue)
caxiscen;
colorbar
title('G - eastward velocity')
ylim([0,12])
hold on
contour(datetime_events*24,12-zsG_isle,squeeze(nanmean(tG_events,3))',[0:.1:100],'color',[0.5 0.5 0.5]), shading flat
[c,h] = contour(datetime_events*24,12-zsG_isle,squeeze(nanmean(tG_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2);
hold off
caxis([-0.15,0.15])
ylabel('[mab]')



subplot(15,2,[18:2:20])
pcolorjw(mitime_events*24,G.M.z(zi),squeeze(nanmean(imag(wg_events(:,zi,:)),3))')
%pcolorjw(mitime_events,G.M.z(zi),squeeze(nanmean(imag(wg_events(:,zi,:)),3))')
shading flat
colormap(redblue)
caxiscen;
colorbar
title('G - northward velocity')
ylim([0,12])
hold on
contour(datetime_events*24,12-zsG_isle,squeeze(nanmean(tG_events,3))',[0:.1:100],'color',[0.5 0.5 0.5]), shading flat
[c,h] = contour(datetime_events*24,12-zsG_isle,squeeze(nanmean(tG_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2);
hold off
caxis([-0.05,0.05])
ylabel('[mab]')


zi = 1:15;
subplot(15,2,[23:2:29])
pcolorjw(mitime_events*24,I.M.z(zi),squeeze(nanmean(real(wi_events(:,zi,:)),3))')
%pcolorjw(mitime_events,I.M.z(zi),squeeze(nanmean(real(wi_events(:,zi,:)),3))')
shading flat
colormap(redblue)
caxiscen;
colorbar
caxis([-0.15,0.15])
ylim([0,22])
hold on 
contour(datetime_events*24,22-zsI_isle,squeeze(nanmean(tI_events,3))',[0:.1:100],'color',[0.5 0.5 0.5]), shading flat
contour(datetime_events*24,22-zsI_isle,squeeze(nanmean(tI_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2), shading flat
hold off
title('I - eastward velocity')
xlabel('time relative to event start [h]')
ylabel('[mab]')


subplot(15,2,[24:2:30])
pcolorjw(mitime_events*24,I.M.z(zi),squeeze(nanmean(imag(wi_events(:,zi,:)),3))')
%pcolorjw(mitime_events,I.M.z(zi),squeeze(nanmean(imag(wi_events(:,zi,:)),3))')
shading flat
colormap(redblue)
caxiscen;
colorbar
ylim([0,22])
hold on 
contour(datetime_events*24,22-zsI_isle,squeeze(nanmean(tI_events,3))',[0:.1:100],'color',[0.5 0.5 0.5]), shading flat
[c,h] = contour(datetime_events*24,22-zsI_isle,squeeze(nanmean(tI_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2), shading flat
hold off
title('I - northward velocity')
caxis([-0.05,0.05])
xlabel('time relative to event start [h]')
ylabel('[mab]')

print -dpng figures_paper/fig_h_event_composite_allsites.png

