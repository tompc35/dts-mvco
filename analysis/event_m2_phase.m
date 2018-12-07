clear all

load_dts_isle_data
detide_h;
m2_freq = freq(5);

m2_cycle = m2_amp(3)*exp(i*2*pi*(m2_freq*24*datetime-m2_pha(3)));
phase = angle(m2_cycle);

% interpolate seagauge temp at C to DTS time base
tsgc = interp1(mdaysg,wtsg(:,3),datetime);

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

% combine temperature mooring and cal temp at site H
tHm = interp1(mday_isle,wtH_isle,datetime);
tH = [tHm,tcal3];
zsH_isle = [zsI_isle, 16];

%%
Huint = interp1(H.ttime,H.uu(2,:),datetime);
Cuint = interp1(C.M.mtime,C.M.evm(12,:),datetime);

nfilt = 3;
thresh = -0.25;
[eventih,event_daten] = get_event_indices_dTdt(boxfilt(tH(:,end),nfilt),datetime,thresh);
[eventic,event_daten] = get_event_indices_dTdt(boxfilt(tsgc,nfilt),datetime,thresh);
[eventii,event_daten] = get_event_indices_dTdt(boxfilt(tI(:,end-3),nfilt),datetime,thresh);
[eventie,event_daten] = get_event_indices_dTdt(boxfilt(te(:,end),nfilt),datetime,thresh);

%%
close all
figure

subplot(211)
plot(H.ttime,H.uu(3,:))
hold on
plot(H.ttime,H.utide(3,:))
hold on
plot(datetime,m2_cycle)
xlim([datenum('8-1-2014'),datenum('9-5-2014')])
hold off

subplot(212)
plot(datetime,angle(m2_cycle));
xlim([datenum('8-1-2014'),datenum('9-5-2014')])

%%
figure
set(gcf, 'PaperSize', [13.0 13.0]);
set(gcf, 'PaperPosition', [0 0 8.5 4.0])

subplot(1,2,1)
he  = rose(phase(eventie),12);
set(he,'Linestyle','none','Color','w','linewidth',1)
hold on
hh  = rose(phase(eventih),12);
set(hh,'Linestyle','-','Color','r','linewidth',3)
hc = rose(phase(eventic),12);
set(hc,'Linestyle','-','Color','b','linewidth',3)
hold off
leg = legend([hh,hc],'site H (bottom)','site C (bottom)')
pos = get(leg,'Position')
set(leg,'fontsize',12,'Position',[pos(1)+0.1 pos(2)+0.066 pos(3) pos(4)])
text(11,3,'eastward','fontsize',14)
text(-26,3,'westward','fontsize',14)

subplot(1,2,2)
hh  = rose(phase(eventie),12);
set(hh,'Linestyle','-','Color','r','linewidth',3)
hold on
hc = rose(phase(eventii),12);
set(hc,'Linestyle','-','Color','b','linewidth',3)
hold off
leg = legend('site E (bottom)','site I (13 mab)')
pos = get(leg,'Position')
set(leg,'fontsize',12,'Position',[pos(1)+0.1 pos(2)+0.066 pos(3) pos(4)])
print('-dpng','../figures/fig_event_polar_histogram')


%% Analyze high pass filtered data
datei = find(datetime >= datenum('13 July 2014') & datetime < datenum('5 September 2014'));

dTdte = 0.5*(te(datei(3):datei(end),end)-te(datei(1):datei(end-2),end))/(600);
dTdth = 0.5*(tcal3(datei(3):datei(end),end)-tcal3(datei(1):datei(end-2),end))/(600);

nanstd(te(datei,end)-pl64tc(te(datei,end)))
nanstd(tcal3(datei,end)-pl64tc(tcal3(datei,end)))
nanstd(tsgc(datei,end)-pl64tc(tsgc(datei,end)))
nanstd(tI(datei,end)-pl64tc(tI(datei,end)))