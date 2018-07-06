% script calibrate_dts_linearsolve_gammathresh_distavg_writenc.m
% --------------------------------------------------------------
% DTS calibration using full (time-averaged) data set.
% 
% Flag points based on a threshold of the "gamma"  calibration parameter,
% since low gamma values are associated with large errors.
%
% Average to 5m, trim anchor and tower segment, and assign lat/lon points.
%
% Write data to a NetCDF file.
%
% Stokes and Anti-Stokes backscatter is best averaged in distance as the log power
% ratio (LPR):
% log(Stokes/AntiStokes) 
%
% Tom Connolly
% January 2015

clear all
run ../data_paths

config_name = 'asit';

out_dir = dts_dir;

dts_nc = [dts_dir 'DTS' config_name '_timeavg.nc'];
cal_nc = [dts_dir 'DTS' config_name '_cal.nc'];

out_nc = [out_dir 'DTS' config_name '_proc.nc']; 

wtcalfile = [cal_dir 'WaterTempPro/temppro_cals.mat'];
gps_file = [cal_dir 'GPS/DTS_GPS_westcable.txt'];

gamma_thresh = 300;     % threshold for gamma flag
dz = 5;                % averaging distance (m)
zs = 160;               % start distance - end of anchor coil
                        % *** this is also where the GPS positions start
                        % ***
                        
%%

distance = ncread(dts_nc,'distance');
datetime = ncread(dts_nc,'datetime');

% add datetime offset if necessary
nci = ncinfo(dts_nc);
if ~isempty(strmatch('datetime_offset',{nci.Variables.Name}))
    datetime_offset = ncread(dts_nc,'datetime_offset');
    datetime = datetime - datetime_offset/24; % adjust all times to GMT
end

% Load water temp pro calibration info
wt = load(wtcalfile);

% Load GPS data
[lon_gps,lat_gps,alt_gps,jnk1,jnk2,waypt,date_gps,time_gps] = ...
    textread(gps_file,'%n%n%n%s%s%s%s%s','headerlines',1,'delimiter',',');
% use only points with time info
di = find(~strcmp(date_gps,''));
lon_gps = lon_gps(di);
lat_gps = lat_gps(di);
date_gps = date_gps(di);
time_gps = time_gps(di);
% calculate distance relative to ASIT (start point of GPS data)
lon_asit = -70.5667;
lat_asit = 41.3250;
[x,y] = latlon2xy([lat_asit; lat_gps],[lon_asit; lon_gps],lat_asit,lon_asit);
dist_gps = cumsum([sqrt(diff(x).^2 + diff(y).^2)])*1000;
% select only unique distances so that vector is monotonic and increasing 
[dist_gps,iu,~] = unique(dist_gps);
lon_gps = lon_gps(iu);
lat_gps = lat_gps(iu);

%%

% specify indices of the reference sections
% ~11m offset between distances from splice and distance detected by
% instrument

refzi1 = 15:58;     % cal bath #1 (cooler)
refzi2 = 65:92;     % cal bath #2 (tote)
refzi3 = find(distance >= 132 & distance <= 154);   % anchor / SBE39 0650
refzi4 = find(distance >= 920 & distance <= 930); % WT Pro - 1269446 
refzi5 = find(distance >= 3878 & distance <= 4004); % WT Pro - 1269445 
refzi6 = find(distance >= 4572 & distance <= 4582); % WT Pro - 1269447 

%%% Load calibration data %%%
tcal1s = ncread(cal_nc,'t_687'); % Seabird
tcal1b = ncread(dts_nc,'tref_2'); % CTEMPS sensor
hi = find(isfinite(tcal1s+tcal1b));
p = polyfit(tcal1b(hi),tcal1s(hi),1);
tcal1_itrp = p(1)*tcal1b + p(2);
tcal1 = tcal1s;
ni = find(isnan(tcal1s));
tcal1(ni) = tcal1_itrp(ni);

tcal2 = ncread(cal_nc,'t_648'); % Seabird
tcal2b = ncread(dts_nc,'tref_1'); % CTEMPS sensor

tcal3 = ncread(cal_nc,'t_650'); % Seabird

% Water Temp Pros (apply calibration)
tcal4p = ncread(cal_nc,'t_446');
j = find(wt.id==1269446);
tcal4=polyval(wt.Tcal(j(end),:),tcal4p);

tcal5p = ncread(cal_nc,'t_445');
j = find(wt.id==1269445);
tcal5=polyval(wt.Tcal(j(end),:),tcal5p);

tcal6p = ncread(cal_nc,'t_447');
j = find(wt.id==1269447);
tcal6=polyval(wt.Tcal(j(end),:),tcal6p);

tint = ncread(dts_nc,'tref_int');

%%% Load DTS data %%%
z1 = mean(ncread(dts_nc,'distance',refzi1(1),length(refzi1)));
stokes1 = ncread(dts_nc,'Stokes',[refzi1(1) 1],[length(refzi1) Inf])';
antistokes1 = ncread(dts_nc,'AntiStokes',[refzi1(1) 1],[length(refzi1) Inf])';
lpr1 = mean(log(stokes1.*antistokes1.^-1),2); % log power ratio, averaged over reference section

z2 = mean(ncread(dts_nc,'distance',refzi2(1),length(refzi2)));
stokes2 = ncread(dts_nc,'Stokes',[refzi2(1) 1],[length(refzi2) Inf])';
antistokes2 = ncread(dts_nc,'AntiStokes',[refzi2(1) 1],[length(refzi2) Inf])';
lpr2 = mean(log(stokes2.*antistokes2.^-1),2); % log power ratio, averaged over reference section

z3 = mean(ncread(dts_nc,'distance',refzi3(1),length(refzi3)));
stokes3 = ncread(dts_nc,'Stokes',[refzi3(1) 1],[length(refzi3) Inf])';
antistokes3 = ncread(dts_nc,'AntiStokes',[refzi3(1) 1],[length(refzi3) Inf])';
lpr3 = mean(log(stokes3.*antistokes3.^-1),2); % log power ratio, averaged over reference section

z4 = mean(ncread(dts_nc,'distance',refzi4(1),length(refzi4)));
stokes4 = ncread(dts_nc,'Stokes',[refzi4(1) 1],[length(refzi4) Inf])';
antistokes4 = ncread(dts_nc,'AntiStokes',[refzi4(1) 1],[length(refzi4) Inf])';
lpr4 = mean(log(stokes4.*antistokes4.^-1),2); % log power ratio, averaged over reference section

z5 = mean(ncread(dts_nc,'distance',refzi5(1),length(refzi5)));
stokes5 = ncread(dts_nc,'Stokes',[refzi5(1) 1],[length(refzi5) Inf])';
antistokes5 = ncread(dts_nc,'AntiStokes',[refzi5(1) 1],[length(refzi5) Inf])';
lpr5 = mean(log(stokes5.*antistokes5.^-1),2); % log power ratio, averaged over reference section

z6 = mean(ncread(dts_nc,'distance',refzi6(1),length(refzi6)));
stokes6 = ncread(dts_nc,'Stokes',[refzi6(1) 1],[length(refzi6) Inf])';
antistokes6 = ncread(dts_nc,'AntiStokes',[refzi6(1) 1],[length(refzi6) Inf])';
lpr6 = mean(log(stokes6.*antistokes6.^-1),2); % log power ratio, averaged over reference section

%%

gamma = nan(size(datetime));
C = nan(size(datetime));
deltaalpha = nan(size(datetime));

disp('calculating calibration parameters')
for ii = 1:length(datetime)
   A = [1  -tcal1(ii)-273.15  (tcal1(ii)+273.15)*z1/1000;
        1  -tcal3(ii)-273.15  (tcal3(ii)+273.15)*z3/1000;
        1  -tcal4(ii)-273.15  (tcal4(ii)+273.15)*z4/1000;
        1  -tcal5(ii)-273.15  (tcal5(ii)+273.15)*z5/1000;
        1  -tcal6(ii)-273.15  (tcal6(ii)+273.15)*z6/1000];  

   b =  [(tcal1(ii)+273.15)*lpr1(ii);
         (tcal3(ii)+273.15)*lpr3(ii);
         (tcal4(ii)+273.15)*lpr4(ii);
         (tcal5(ii)+273.15)*lpr5(ii);
         (tcal6(ii)+273.15)*lpr6(ii)];
     
   gi = find(isfinite(b));
   b = b(gi,:);
   A = A(gi,:);
     
   % solve A\b and calculate determinant of A matrix
   if length(b) < 3
        x = [NaN; NaN; NaN];
   else
        x = A\b;
   end
   
   gamma(ii) = x(1);
   C(ii) = x(2);
   deltaalpha(ii) = x(3);
   
end

% time indices with low gamma values
ig = find(abs(gamma)<gamma_thresh);

%%

% Define distance variable for calibrated, distance-averaged temperature
zbin = zs+dz/2:dz:distance(end);
z = (zbin(1:end-1)+zbin(2:end))/2;

lon = interp1(dist_gps,lon_gps,z-zs);
lat = interp1(dist_gps,lat_gps,z-zs);

% % Create a new netcdf file for calibrated, distance-averaged temperature
% ncid = netcdf.create(out_nc,'CLOBBER');

delete(out_nc)

nccreate(out_nc,'tempC','Dimensions',{'x',length(z),'t',Inf});
ncwriteatt(out_nc,'tempC','variable','Temperature along cable (calibrated)')
ncwriteatt(out_nc,'tempC','units','degrees Celcius')

nccreate(out_nc,'lpr','Dimensions',{'x','t'});
ncwriteatt(out_nc,'lpr','variable','log power ratio, ln(Stokes/AntiStokes)')
ncwriteatt(out_nc,'lpr','units','none')

nccreate(out_nc,'distance','Dimensions',{'x'});
ncwriteatt(out_nc,'distance','variable','distance (relative to instrument)')
ncwriteatt(out_nc,'distance','units','meters')
ncwrite(out_nc,'distance',z)

nccreate(out_nc,'lon','Dimensions',{'x'});
ncwriteatt(out_nc,'lon','variable','longitude')
ncwriteatt(out_nc,'lon','units','decimal degrees east')
ncwrite(out_nc,'lon',lon)

nccreate(out_nc,'lat','Dimensions',{'x'});
ncwriteatt(out_nc,'lat','variable','longitude')
ncwriteatt(out_nc,'lat','units','decimal degrees east')
ncwrite(out_nc,'lat',lat)

nccreate(out_nc,'datetime','Dimensions',{'t'});
ncwriteatt(out_nc,'datetime','variable','Serial Date Number, GMT (Matlab convention)')
ncwriteatt(out_nc,'datetime','units','days')
ncwriteatt(out_nc,'datetime','note','A serial date number of 1 corresponds to Jan-1-0000')
ncwrite(out_nc,'datetime',datetime)

nccreate(out_nc,'gamma','Dimensions',{'t'});
ncwriteatt(out_nc,'gamma','variable','calibration parameter, gamma')
ncwriteatt(out_nc,'gamma','units','degrees Kelvin')
ncwrite(out_nc,'gamma',gamma)

nccreate(out_nc,'C','Dimensions',{'t'});
ncwriteatt(out_nc,'C','variable','calibration parameter, C')
ncwriteatt(out_nc,'C','units','none')
ncwrite(out_nc,'C',C)

nccreate(out_nc,'deltaalpha','Dimensions',{'t'});
ncwriteatt(out_nc,'deltaalpha','variable','calibration parameter, deltaalpha')
ncwriteatt(out_nc,'deltaalpha','units','m^-1')
ncwrite(out_nc,'deltaalpha',deltaalpha)

% Loop through distance  segments, calculate mean log power ratio
% and use calibration coefficients to calculate temperature


disp('writing new file')
for ii = 1:[length(z)]
    if mod(ii,10)==0
        disp(['iteration ' num2str(ii) '/' num2str(length(z)-1)])
    end
    
    % Find distance indices of this section
    zi = find(distance >= zbin(ii) & distance < zbin(ii+1));
    
    % Load DTS data from this section
    z_sec = mean(ncread(dts_nc,'distance',zi(1),length(zi)));
    stokes_sec = ncread(dts_nc,'Stokes',[zi(1) 1],[length(zi) Inf])';
    antistokes_sec = ncread(dts_nc,'AntiStokes',[zi(1) 1],[length(zi) Inf])';
    lpr_sec = [mean(log(stokes_sec.*antistokes_sec.^-1),2)]'; % log power ratio, averaged over reference section
    
    % calculate temperature 
    t_sec = gamma'.*(lpr_sec + C' - deltaalpha'*z_sec/1000).^-1 - 273.15;
    
    % apply flags
    t_sec(ig) = NaN;
        
    ncwrite(out_nc,'lpr',lpr_sec,[ii 1],[1 1]);
    ncwrite(out_nc,'tempC',t_sec,[ii 1],[1 1]);
end

ncwriteatt(out_nc,'/','source',['Created by calibrate_dts_linearsolve_gammathresh_distavg_writenc_cal.m, ' datestr(now)])
ncwriteatt(out_nc,'/','calibration pts',['all available'])

