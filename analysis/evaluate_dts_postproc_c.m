clear all
close all

plot_event = 1;

load_dts_isle_data
nfilt = 3;

% interpolate seagauge temp at C to DTS time base
tsgc = interp1(mdaysg,wtsg(:,3),datetime);
[eventi,event_daten] = get_event_indices_dTdt(boxfilt(tsgc,nfilt),datetime,-0.5);

tsgcf = nan(size(tsgc));
tdtsf = nan(size(tsgc));
tsgcf(2:end-1) = boxfilt(tsgc,nfilt);
tdtsf(2:end-1) = boxfilt(tempC(zic,:),nfilt);

% tsgcf = boxfilt(tsgc,nfilt);
% tdtsf = boxfilt(tempC(zic,:),nfilt);

% Full deployment
rmse = sqrt(nanmean((tsgc-tempC(zic,:)').^2));
bias = nanmean(tsgc-tempC(zic,:)');
disp('C bottom temperature comparison')
disp('-------------------------------')
disp('Entire deployment:')
disp(['RMSE = ' num2str(rmse,2)])
disp(['bias = ' num2str(bias,2)])

% Stratified season
ti = find(datetime<datenum('Sept-6-2014'));
rmse = sqrt(nanmean((tsgc(ti)-tempC(zic,ti)').^2));
bias = nanmean(tsgc(ti)-tempC(zic,ti)');
disp('Stratified season (through 6 Sep 2014):')
disp(['RMSE = ' num2str(rmse,2)])
disp(['bias = ' num2str(bias,2)])

% Unstratified season
ti = find(datetime>=datenum('Sept-6-2014'));
rmse = sqrt(nanmean((tsgc(ti)-tempC(zic,ti)').^2));
bias = nanmean(tsgc(ti)-tempC(zic,ti)');
disp('Unstratified season (after 6 Sep 2014):')
disp(['RMSE = ' num2str(rmse,2)])
disp(['bias = ' num2str(bias,2)])

% distance between points
[x,y] = latlon2xy(lat_isle(3),lon_isle(3),lat_dts(zic),lon_dts(zic));
d = sqrt(x^2+y^2);
disp(['distance between points = ' num2str(d*1000,3) 'm'])
disp(' ')

% dT/dt comparison
tc = tempC(zic,:)';
dt = (datetime(2)-datetime(1))*24;
dTdtC = nan(size(tsgcf));
dTdtC(2:end-1) = (tdtsf(3:end)-tdtsf(1:end-2))/(2*dt);
dTdtcalC = nan(size(tsgcf));
dTdtcalC(2:end-1) = (tsgcf(3:end)-tsgcf(1:end-2))/(2*dt);

%%
close all

figure
plot(datetime,tsgcf)
hold on, plot(datetime,tdtsf)
datetick('x')

figure
plot(tsgc,tempC(zic,:),'.')

figure
plot(dTdtC,dTdtcalC,'.')
gi = find(isfinite(dTdtC+dTdtcalC));
[r,p] = corrcoef(dTdtC(gi),dTdtcalC(gi));

disp(['dT/dt r = ' num2str(r(2))])