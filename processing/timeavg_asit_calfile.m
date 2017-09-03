% timeavg_asit_dts_ncfile

nc_dir = '/media/tompc/data/DTS_nc/';
nc_out = [nc_dir 'DTSasit_cal.nc'];

t_start = datenum('11-Jul-2014 18:20:00');
t_end = datenum('31-Oct-2014 13:30:00');
t_interval = 10*60/86400; % time interval (days)
datetime_new = t_start:t_interval:t_end;

%%% load variables from seabird sensors

sbe_dir = '/media/tompc/data/DTS_cal/Seabird/';

sbe_name = 'DTS_temp_sbe0648_11112014.asc';
[tempC_648,day,time] = textread([sbe_dir sbe_name],'%n%s%s','headerlines',36,'delimiter',',');
dnum_648 = datenum([char(day) char(time)],'dd mmm yyyyHH:MM:SS');

sbe_name = 'DTS_temp_sbe0650_11112014.asc';
[tempC_650,day,time] = textread([sbe_dir sbe_name],'%n%s%s','headerlines',36,'delimiter',',');
dnum_650 = datenum([char(day) char(time)],'dd mmm yyyyHH:MM:SS');

sbe_name = 'DTS_temp_sbe0652_11112014.asc';
[tempC_652,day,time] = textread([sbe_dir sbe_name],'%n%s%s','headerlines',36,'delimiter',',');
dnum_652 = datenum([char(day) char(time)],'dd mmm yyyyHH:MM:SS');

sbe_name = 'DTS_temp_sbe0687_11112014.asc';
[tempC_687,day,time] = textread([sbe_dir sbe_name],'%n%s%s','headerlines',36,'delimiter',',');
dnum_687 = datenum([char(day) char(time)],'dd mmm yyyyHH:MM:SS');

%%% load variables from water temp pro sensors

wtpro_dir = '/media/tompc/data/DTS_cal/WaterTempPro/';

wtpro_name = 'DTS_cal_1038035_20141113.csv';
[num,daytime,tempC_035] = textread([wtpro_dir wtpro_name],'%n%s%n','headerlines',2,'delimiter',',');
dnum_035 = datenum([char(daytime)],'mm/dd/yyyy HH:MM:SS');

wtpro_name = 'DTS_cal_1269445_20141113.csv';
[num,daytime,tempC_445] = textread([wtpro_dir wtpro_name],'%n%s%n','headerlines',2,'delimiter',',');
dnum_445 = datenum([char(daytime)],'mm/dd/yyyy HH:MM:SS');

wtpro_name = 'DTS_cal_1269446_20141113.csv';
[num,daytime,tempC_446] = textread([wtpro_dir wtpro_name],'%n%s%n','headerlines',2,'delimiter',',');
dnum_446 = datenum([char(daytime)],'mm/dd/yyyy HH:MM:SS');

wtpro_name = 'DTS_cal_1269447_20141113.csv';
[num,daytime,tempC_447] = textread([wtpro_dir wtpro_name],'%n%s%n','headerlines',2,'delimiter',',');
dnum_447 = datenum([char(daytime)],'mm/dd/yyyy HH:MM:SS');

%%% Set up NetCDF file

% create new file for output, overwriting existing file
ncid = netcdf.create(nc_out,'CLOBBER');

% define time dimension with unlimited length (t)
t_dimID = netcdf.defDim(ncid,'t',...
    netcdf.getConstant('NC_UNLIMITED'));

% define variables
date_varID = netcdf.defVar(ncid,'datetime','NC_DOUBLE',t_dimID);
netcdf.putAtt(ncid,date_varID,'variable','Serial Date Number, GMT (Matlab convention)')
netcdf.putAtt(ncid,date_varID,'units','days')
netcdf.putAtt(ncid,date_varID,'note','A serial date number of 1 corresponds to Jan-1-0000')

tint_varID = netcdf.defVar(ncid,'t_648','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  
netcdf.putAtt(ncid,tint_varID,'instrument','Seabird 39 S/N 0648')  
netcdf.putAtt(ncid,tint_varID,'location','tower calibration bath - blue tote')  

tint_varID = netcdf.defVar(ncid,'t_650','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius') 
netcdf.putAtt(ncid,tint_varID,'instrument','Seabird 39 S/N 0650')  
netcdf.putAtt(ncid,tint_varID,'location','west cable - tower leg')

tint_varID = netcdf.defVar(ncid,'t_652','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  
netcdf.putAtt(ncid,tint_varID,'instrument','Seabird 39 S/N 0652')  
netcdf.putAtt(ncid,tint_varID,'location','east cable - tower leg')

tint_varID = netcdf.defVar(ncid,'t_687','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  
netcdf.putAtt(ncid,tint_varID,'instrument','Seabird 39 S/N 0687')  
netcdf.putAtt(ncid,tint_varID,'location','tower calibration bath - cooler')

tint_varID = netcdf.defVar(ncid,'t_446','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  
netcdf.putAtt(ncid,tint_varID,'instrument','Water Temp Pro S/N 1269446')  
netcdf.putAtt(ncid,tint_varID,'location','west cable')

tint_varID = netcdf.defVar(ncid,'t_445','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  
netcdf.putAtt(ncid,tint_varID,'instrument','Water Temp Pro S/N 1269445')  
netcdf.putAtt(ncid,tint_varID,'location','west cable')

tint_varID = netcdf.defVar(ncid,'t_447','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  
netcdf.putAtt(ncid,tint_varID,'instrument','Water Temp Pro S/N 1269447')  
netcdf.putAtt(ncid,tint_varID,'location','west cable')

tint_varID = netcdf.defVar(ncid,'t_035','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  
netcdf.putAtt(ncid,tint_varID,'instrument','Water Temp Pro S/N 1038035')  
netcdf.putAtt(ncid,tint_varID,'location','west cable')

netcdf.endDef(ncid);

nt = length(datetime_new);
netcdf.putVar(ncid,date_varID,0,nt,datetime_new);

netcdf.close(ncid);

sensorid = {'648','650','652','687','446','445','447','035'};

tic 
disp('averaging and writing new file')
for ii = 1:length(datetime_new)
    if mod(ii,100)==0
        disp(['iteration ' num2str(ii) ', ' datestr(datetime_new(ii))])
    end
    
    % Average and write data from each sensor
    for jj = 1:length(sensorid)
        eval(['dnum_tmp = dnum_' char(sensorid(jj)) ';']);
        ti = find(dnum_tmp >= datetime_new(ii) - t_interval/2 & ...
            dnum_tmp < datetime_new(ii) + t_interval/2);
        if isempty(ti)
            t_tmp = NaN;
        else
            eval(['t_tmp = mean(tempC_' char(sensorid(jj)) '(ti));']);
        end
        ncwrite(nc_out,['t_' char(sensorid(jj))],t_tmp,ii,1);
    end
end
toc
