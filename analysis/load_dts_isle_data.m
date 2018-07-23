run ../data_paths

load(microcat_mat)
mday_isle = jd2matday(jd_isle);

load(seagauge_mat)
mdaysg = jd2matday(jdsg');

Equad = load(quadpod_e_mat);

m_nc = [adcp_nc_dir 'ISLE_station_C_velocity_wd14_23-Feb-2017.nc'];
C.M.z = ncread(m_nc,'bin_height');
C.M.mtime = ncread(m_nc,'datetime')+datenum(1970,1,1);
C.M.evm = ncread(m_nc,'East_vel');
C.M.nvm = ncread(m_nc,'North_vel');
C.M.evm(find(C.M.evm == 999)) = NaN;
C.M.nvm(find(C.M.nvm == 999)) = NaN;

m_nc = [adcp_nc_dir 'ISLE_station_E_velocity_wd15_23-Feb-2017.nc'];
E.M.z = ncread(m_nc,'bin_height');
E.M.mtime = ncread(m_nc,'datetime');
E.M.evm = ncread(m_nc,'East_vel');
E.M.nvm = ncread(m_nc,'North_vel');
E.M.evm(find(E.M.evm == 999)) = NaN;
E.M.nvm(find(E.M.nvm == 999)) = NaN;

m_nc = [adcp_nc_dir 'ISLE_station_G_velocity_wd12_23-Feb-2017.nc'];
G.M.z = ncread(m_nc,'bin_height');
G.M.mtime = ncread(m_nc,'datetime')+datenum(1970,1,1);
G.M.evm = ncread(m_nc,'East_vel');
G.M.nvm = ncread(m_nc,'North_vel');
G.M.evm(find(G.M.evm == 999)) = NaN;
G.M.nvm(find(G.M.nvm == 999)) = NaN;

m_nc = [adcp_nc_dir 'ISLE_station_I_velocity_wd22_23-Feb-2017.nc'];
I.M.z = ncread(m_nc,'bin_height');
I.M.mtime = ncread(m_nc,'datetime')+datenum(1970,1,1);
I.M.evm = ncread(m_nc,'East_vel');
I.M.nvm = ncread(m_nc,'North_vel');
I.M.temp = ncread(m_nc,'Temperature');
I.M.evm(find(I.M.evm == 999)) = NaN;
I.M.nvm(find(I.M.nvm == 999)) = NaN;

m_nc = [adcp_nc_dir 'ISLE_station_H_velocity_wd16_23-Feb-2017.nc'];
H.zz = ncread(m_nc,'bin_height');
H.ttime = ncread(m_nc,'datetime')+datenum(1970,1,1);
H.uu = ncread(m_nc,'East_vel');
H.vv = ncread(m_nc,'North_vel');
H.temp = ncread(m_nc,'Temperature');
H.uu(find(H.uu == 999)) = NaN;
H.vv(find(H.vv == 999)) = NaN;

distance = ncread(dts_nc,'distance');
datetime = ncread(dts_nc,'datetime');
tempC = ncread(dts_nc,'tempC');
[~,zic] = min(abs(distance - 4350)); % corner

tcal3 = ncread(cal_nc,'t_650'); % Seabird

% Water Temp Pros (apply calibration)
wt = load(wtcalfile);

tcal4p = ncread(cal_nc,'t_446');
j = find(wt.id==1269446);
tcal4=polyval(wt.Tcal(j(end),:),tcal4p);

tcal5p = ncread(cal_nc,'t_445');
j = find(wt.id==1269445);
tcal5=polyval(wt.Tcal(j(end),:),tcal5p);

tcal6p = ncread(cal_nc,'t_447');
j = find(wt.id==1269447);
tcal6=polyval(wt.Tcal(j(end),:),tcal6p);

clear tcal*p j

% DTS corner index
[~,zic] = min(abs(distance - 4350)); 

% DTS 
lon_dts = ncread(dts_nc,'lon');
lat_dts = ncread(dts_nc,'lat');

% positions of calibration points
d_dts = ncread(dts_nc,'distance');

% calibration points 
% (distances from calibrate_dts_linearsolve_gammathresh_towercal_cal6.m)
d3 = 143.0352;
d4 = 924.9522;
d5 = 3.9433e+03;
d6 = 4.5747e+03;

zi = find(~isnan(lon_dts));

lon3 = lon_dts(min(zi)); % First point is in coil before lat/lon points start
lat3 = lat_dts(min(zi));
lon4 = interp1(d_dts,lon_dts,d4);
lat4 = interp1(d_dts,lat_dts,d4);
lon5 = interp1(d_dts,lon_dts,d5);
lat5 = interp1(d_dts,lat_dts,d5);
lon6 = interp1(d_dts,lon_dts,d6);
lat6 = interp1(d_dts,lat_dts,d6);

clear zi