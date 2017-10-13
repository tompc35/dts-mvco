clear all

load_dts_isle_data

em_dir = '/Users/tomconnolly/work/Data/ISLE/archive_Isle_mooring/'
em_nc = [em_dir 'ISLE_station_E_velocity_wd15_23-Feb-2017.nc'];
E.M.z = ncread(em_nc,'bin_height')
E.M.mtime = ncread(em_nc,'datetime')
E.M.evm = ncread(em_nc,'East_vel')
E.M.nvm = ncread(em_nc,'North_vel')
E.M.evm(find(E.M.evm == 999)) = NaN;
E.M.nvm(find(E.M.nvm == 999)) = NaN;

%%

nfilt = 3;
[eventi,event_daten] = get_event_indices_dTdt(boxfilt(tcal3,nfilt),datetime,-0.5);

%%

% find where H time series starts
ti = 1:length(eventi);

th = interp1(mday_isle,wtH_isle,datetime);
tc = interp1(mday_isle,wtC_isle,datetime);
tG = interp1(mday_isle,wtG_isle,datetime);

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

        tdts_events = nan([length(di) size(tempC,1) length(ti)]);
        wh_events = nan([length(hi) length(H.zz) length(ti)])+i*nan([length(hi) length(H.zz) length(ti)]);
        wi_events = nan([length(Ii) length(I.M.z) length(ti)]);
        wg_events = nan([length(Ii) length(G.M.z) length(ti)]);
        we_events = nan([length(Ii) length(E.M.z) length(ti)]);
        datetime_events = datetime(di)-event_daten(jj);
        mday_isle_events = mday_isle(tim)-event_daten(jj);
        ttime_events = H.ttime(hi)-event_daten(jj);
        mitime_events = I.M.mtime(Ii)-event_daten(jj);
    end
    
    t3_events(:,ii) = tcal3(di);
    th_events(:,:,ii) = th(di,:);
    tc_events(:,:,ii) = tc(di,:);
    te_events(:,:,ii) = te(di,:);
    tI_events(:,:,ii) = tI(di,:);
    tG_events(:,:,ii) = tG(di,:);
    tdts_events(:,:,ii) = tempC(:,di)';
    wh_events(:,:,ii) = H.uu(:,hi)' + i*H.vv(:,hi)';
    wi_events(:,:,ii) = I.M.evm(:,Ii)' + i*I.M.nvm(:,Ii)'; 
    wg_events(:,:,ii) = G.M.evm(:,Ii)' + i*G.M.nvm(:,Ii)'; 
    we_events(:,:,ii) = E.M.evm(:,Ii)' + i*E.M.nvm(:,Ii)'; 
end

t3_events_anom = t3_events-repmat(t3_events(round(length(di)/2),:),[length(di) 1]);

%%
% find events where velocity data exists at H
hidx = find(isfinite(squeeze(sum(wh_events(15,:,:),2))));

% combine temperatures at H
t3_events_tmp(1,:,:) = t3_events;
t_events = [th_events permute(t3_events_tmp,[2 1 3])];
zt_events = [zsH_isle 15];

%%

figure
set(gcf, 'PaperSize', [7.0 7.0]);
set(gcf, 'PaperPosition', [0 0 6.8 6.2])

ax = subplot(15,2,[1:2:5])
pos = get(ax, 'Position');
pos(1) = pos(1)+0.1;
zi = 1:20;
set(ax, 'Position',pos);
pcolorjw(mitime_events*24,E.M.z(zi)+3,squeeze(nanmean(real(we_events(:,zi,:)),3))')
hold on
contour(datetime_events*24,15-zsE_isle,squeeze(nanmean(te_events,3))',[0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events*24,15-zsE_isle,squeeze(nanmean(te_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
title('E eastward velocity')
cbar = colorbar;
set(cbar,'visible','off')
ylabel('[mab]')
ylim([0,15])
caxis([-0.2,0.2])
set(gca,'xtick',[-6:3:6])
xl = xlim;
yl = ylim;
text(xl(1)-0.2*diff(xl),yl(2)+0.15*diff(yl),'a)')

subplot(15,2,[2:2:6])
zi = 1:20;
pcolorjw(mitime_events*24,E.M.z(zi)+3,squeeze(nanmean(imag(we_events(:,zi,:)),3))')
hold on
contour(datetime_events*24,15-zsE_isle,squeeze(nanmean(te_events,3))',[0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events*24,15-zsE_isle,squeeze(nanmean(te_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
title('E northward velocity')
cbar = colorbar;
%set(cbar,'visible','off')
ylabel('[mab]')
ylim([0,15])
caxis([-0.2,0.2])
set(gca,'xtick',[-6:3:6])
xl = xlim;
yl = ylim;
text(xl(1)-0.2*diff(xl),yl(2)+0.15*diff(yl),'b)')

% Standard errors
t_anomstd_e = squeeze(nanstd(te_events(:,end,:)-repmat(nanmean(te_events(:,end,:),1),[length(datetime_events),1,1]),[],3))';
t_anomstderrmax_e = max(t_anomstd_e/(length(isfinite(te_events(1,end,:))))) % divide by number of events and take maximum

u_std_e = squeeze(nanstd(real(we_events(:,1,:)),[],3))';
u_stderrmax_e = max(u_std_e/(length(isfinite(we_events(1,1,:))))) % divide by number of events and take maximum

v_std_e = squeeze(nanstd(imag(we_events(:,1,:)),[],3))';
v_stderrmax_e = max(v_std_e/(length(isfinite(we_events(1,1,:))))) % divide by number of events and take maximum


ax = subplot(15,2,[9:2:13])
pos = get(ax, 'Position');
pos(1) = pos(1)+0.1;
set(ax,'Position',pos)
zi = 1:20;
pcolorjw(ttime_events*24,H.zz,squeeze(nanmean(real(wh_events(:,:,hidx)),3))')
shading flat
colormap(redblue)
caxiscen;
hold on
contour(datetime_events*24,15-zt_events,squeeze(nanmean(t_events(:,:,hidx),3))', ...
            [0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events*24,15-zt_events,squeeze(nanmean(t_events(:,:,hidx),3))', ...
            [0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
cbar1 = colorbar;
cpos1 = get(cbar1,'position');
set(cbar1,'visible','off')
title('H - eastward velocity')
ylim([0,15])
caxis([-0.2,0.2])
ylabel('[mab]')
set(gca,'xtick',[-6:3:6])
xl = double(xlim);
yl = double(ylim);
text(xl(1)-0.2*diff(xl),yl(2)+0.15*diff(yl),'c)')

subplot(15,2,[10:2:14])
pcolorjw(ttime_events*24,H.zz,squeeze(nanmean(imag(wh_events(:,:,hidx)),3))')
%pcolorjw(ttime_events,H.zz,squeeze(nanmean(imag(wh_events(:,:,hidx)),3))')
shading flat
colormap(redblue)
caxiscen;
cbar2 = colorbar;
cpos2 = get(cbar2,'position')
set(cbar2,'visible','off')
hold on
contour(datetime_events*24,15-zt_events,squeeze(nanmean(t_events(:,:,hidx),3))', ...
            [0:.1:100],'color',[0.5 0.5 0.5])
contour(datetime_events*24,15-zt_events,squeeze(nanmean(t_events(:,:,hidx),3))', ...
            [0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2)
title('H - northward velocity')
ylim([0,15])
caxis([-0.2,0.2])
ylabel('[mab]')
set(gca,'xtick',[-6:3:6])
xl = double(xlim);
yl = double(ylim);
text(xl(1)-0.2*diff(xl),yl(2)+0.15*diff(yl),'d)')

% Standard errors
t_anomstd = squeeze(nanstd(t_events(:,end,hidx)-repmat(nanmean(t_events(:,end,hidx),1),[length(datetime_events),1,1]),[],3))';
t_anomstderrmax = max(t_anomstd/(length(isfinite(t_events(1,end,hidx))))) % divide by number of events and take maximum

u_std_h = squeeze(nanstd(real(wh_events(:,end,:)),[],3))';
u_stderrmax_h = max(u_std_h/(length(isfinite(wh_events(1,1,:))))) % divide by number of events and take maximum

v_std_h = squeeze(nanstd(imag(wh_events(:,end,:)),[],3))';
v_stderrmax_h = max(v_std_h/(length(isfinite(we_events(1,1,:))))) % divide by number of events and take maximum


%%% G velocity %%%

zi = 1:15;
ax = subplot(15,2,[17:2:19])
pos = get(ax, 'Position');
pos(1) = pos(1)+0.1;
set(ax,'Position',pos)
pcolorjw(mitime_events*24,G.M.z(zi),squeeze(nanmean(real(wg_events(:,zi,:)),3))')
shading flat
colormap(redblue)
caxiscen;
cbar = colorbar;
set(cbar,'visible','off')
ylim([0,12])
hold on
contour(datetime_events*24,12-zsG_isle,squeeze(nanmean(tG_events,3))',[0:.1:100],'color',[0.5 0.5 0.5]), shading flat
[c,h] = contour(datetime_events*24,12-zsG_isle,squeeze(nanmean(tG_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2);
hold off
caxis([-0.2,0.2])
ylabel('[mab]')
set(gca,'xtick',[-6:3:6])
xl = double(xlim);
yl = double(ylim);
text(xl(1)-0.2*diff(xl),yl(2)+0.15*diff(yl),'e)')

% Standard errors
t_anomstd_g = squeeze(nanstd(tG_events(:,end,:)-repmat(nanmean(tG_events(:,end,:),1),[length(datetime_events),1,1]),[],3))';
t_anomstderrmax_g = max(t_anomstd_g/(length(isfinite(tG_events(1,end,:))))) % divide by number of events and take maximum

u_std_g = squeeze(nanstd(real(wg_events(:,1,:)),[],3))';
u_stderrmax_g = max(u_std_g/(length(isfinite(wg_events(1,1,:))))) % divide by number of events and take maximum

v_std_g = squeeze(nanstd(imag(wg_events(:,1,:)),[],3))';
v_stderrmax_g = max(v_std_g/(length(isfinite(wg_events(1,1,:))))) % divide by number of events and take maximum

subplot(15,2,[18:2:20])
pcolorjw(mitime_events*24,G.M.z(zi),squeeze(nanmean(imag(wg_events(:,zi,:)),3))')
shading flat
colormap(redblue)
caxiscen;
cbar = colorbar;
set(cbar,'visible','off')
title('G - northward velocity')
ylim([0,12])
hold on
contour(datetime_events*24,12-zsG_isle,squeeze(nanmean(tG_events,3))',[0:.1:100],'color',[0.5 0.5 0.5]), shading flat
[c,h] = contour(datetime_events*24,12-zsG_isle,squeeze(nanmean(tG_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2);
hold off
caxis([-0.2,0.2])
ylabel('[mab]')
set(gca,'xtick',[-6:3:6])
xl = double(xlim);
yl = double(ylim);
text(xl(1)-0.2*diff(xl),yl(2)+0.15*diff(yl),'f)')

zi = 1:15;
ax = subplot(15,2,[23:2:29])
pos = get(ax, 'Position');
pos(1) = pos(1)+0.1;
set(ax,'Position',pos)
pcolorjw(mitime_events*24,I.M.z(zi),squeeze(nanmean(real(wi_events(:,zi,:)),3))')
shading flat
colormap(redblue)
caxiscen;
cbar = colorbar;
set(cbar,'visible','off')
caxis([-0.2,0.2])
ylim([0,22])
hold on 
contour(datetime_events*24,22-zsI_isle,squeeze(nanmean(tI_events,3))',[0:.1:100],'color',[0.5 0.5 0.5]), shading flat
contour(datetime_events*24,22-zsI_isle,squeeze(nanmean(tI_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2), shading flat
hold off
title('I - eastward velocity')
xlabel('time relative to event start [h]')
ylabel('[mab]')
set(gca,'xtick',[-6:3:6])
xl = double(xlim);
yl = double(ylim);
text(xl(1)-0.2*diff(xl),yl(2)+0.15*diff(yl),'g)')

subplot(15,2,[24:2:30])
pcolorjw(mitime_events*24,I.M.z(zi),squeeze(nanmean(imag(wi_events(:,zi,:)),3))')
shading flat
colormap(redblue)
caxiscen;
cbar = colorbar;
set(cbar,'visible','off')
ylim([0,22])
hold on 
contour(datetime_events*24,22-zsI_isle,squeeze(nanmean(tI_events,3))',[0:.1:100],'color',[0.5 0.5 0.5]); shading flat
[c,h] = contour(datetime_events*24,22-zsI_isle,squeeze(nanmean(tI_events,3))',[0:.5:100],'color',[0.3 0.3 0.3],'linewidth',2); shading flat
hold off
title('I - northward velocity')
caxis([-0.2,0.2])
xlabel('time relative to event start [h]')
ylabel('[mab]')
set(gca,'xtick',[-6:3:6]);
xl = double(xlim);
yl = double(ylim);
text(xl(1)-0.2*diff(xl),yl(2)+0.15*diff(yl),'h)')

% Standard errors
t_anomstd_i = squeeze(nanstd(tI_events(:,end,:)-repmat(nanmean(tI_events(:,end,:),1),[length(datetime_events),1,1]),[],3))';
t_anomstderrmax_i = max(t_anomstd_i/(length(isfinite(tI_events(1,end,:))))) % divide by number of events and take maximum

u_std_i = squeeze(nanstd(real(wi_events(:,1,:)),[],3))';
u_stderrmax_i = max(u_std_i/(length(isfinite(wi_events(1,1,:))))) % divide by number of events and take maximum

v_std_i = squeeze(nanstd(imag(wi_events(:,1,:)),[],3))';
v_stderrmax_i = max(v_std_i/(length(isfinite(wi_events(1,1,:))))) % divide by number of events and take maximum

print -dpng ../figures/fig_h_event_composite_allsites.png