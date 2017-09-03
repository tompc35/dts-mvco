load_data_paths

load(microcat_mat)
mday_isle = jd2matday(jd_isle);

load(seagauge_mat)
mdaysg = jd2matday(jdsg');

load(node_mat)
Equad = load(quadpod_e_mat);

C = load(adcp_c_mat);
G = load(adcp_g_mat);
H = load(adcp_h_mat);
I = load(adcp_i_mat);

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