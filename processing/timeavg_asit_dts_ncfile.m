% timeavg_asit_dts_ncfile

nc_dir = '/media/tompc/data/DTS_nc/';

dts_filename_in = 'DTSasit_01_channel1';

if strcmp(dts_filename_in,'DTSasit_02c_channel1')
    t_start = datenum('17-Jul-2014 19:50:00');
    t_end = datenum('31-Oct-2014 13:30:00');
elseif strcmp(dts_filename_in,'DTSasit_01_channel1')
    t_start = datenum('11-Jul-2014 18:10:00');
    t_end = datenum('17-Jul-2014 19:20:00');    
end

t_interval = 10*60/86400; % time interval (days)

%%%

dts_nc_in = [nc_dir dts_filename_in '.nc'];
dts_nc_out = [nc_dir dts_filename_in '_timeavg.nc'];

datetime = ncread(dts_nc_in,'datetime');
distance = ncread(dts_nc_in,'distance');

datetime_new = t_start:t_interval:t_end;

% create new file for output, overwriting existing file
ncid = netcdf.create(dts_nc_out,'CLOBBER');

% define distance dimension (x)
x_dimID = netcdf.defDim(ncid,'x',length(distance));

% define time dimension with unlimited length (t)
t_dimID = netcdf.defDim(ncid,'t',...
    netcdf.getConstant('NC_UNLIMITED'));

% define variables
dist_varID = netcdf.defVar(ncid,'distance','NC_DOUBLE',x_dimID);
netcdf.putAtt(ncid,dist_varID,'variable','along-cable distance')
netcdf.putAtt(ncid,dist_varID,'units','meters')

date_varID = netcdf.defVar(ncid,'datetime','NC_DOUBLE',t_dimID);
netcdf.putAtt(ncid,date_varID,'variable','Serial Date Number, GMT (Matlab convention)')
netcdf.putAtt(ncid,date_varID,'units','days')
netcdf.putAtt(ncid,date_varID,'note','A serial date number of 1 corresponds to Jan-1-0000')

AntiStokes_varID = netcdf.defVar(ncid,'AntiStokes','NC_DOUBLE',[x_dimID, t_dimID]);
netcdf.putAtt(ncid,AntiStokes_varID,'variable','anti-Stokes intensity')
netcdf.putAtt(ncid,AntiStokes_varID,'units','arbitrary')

Stokes_varID = netcdf.defVar(ncid,'Stokes','NC_DOUBLE',[x_dimID, t_dimID]);
netcdf.putAtt(ncid,Stokes_varID,'variable','Stokes intensity')
netcdf.putAtt(ncid,Stokes_varID,'units','arbitrary')

tint_varID = netcdf.defVar(ncid,'tref_int','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,tint_varID,'variable','Internal Instrument Temperature')
netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  

t1_varID = netcdf.defVar(ncid,'tref_1','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,t1_varID,'variable','Reference Temperature #1')
netcdf.putAtt(ncid,t1_varID,'units','degrees Celcius')

t2_varID = netcdf.defVar(ncid,'tref_2','NC_DOUBLE',[t_dimID]);
netcdf.putAtt(ncid,t2_varID,'variable','Reference Temperature #2')
netcdf.putAtt(ncid,t2_varID,'units','degrees Celcius')   

netcdf.endDef(ncid);

nx = length(distance);
nt = length(datetime_new);
netcdf.putVar(ncid,dist_varID,0,nx,distance);
netcdf.putVar(ncid,date_varID,0,nt,datetime_new);

netcdf.close(ncid);

tic 
disp('averaging and writing new file')
for ii = 1:length(datetime_new)
    if mod(ii,100)==0
        disp(['iteration ' num2str(ii) ', ' datestr(datetime_new(ii))])
    end
    ti = find(datetime >= datetime_new(ii) - t_interval/2 & ...
        datetime < datetime_new(ii) + t_interval/2);
    
    % Average data where available or fill with NaNs
    if isempty(ti)
        Stokes_tmp = nan(nx,1);
        AntiStokes_tmp = nan(nx,1);
        tref_int_tmp = NaN;
        tref_1_tmp = NaN;
        tref_2_tmp = NaN;
    else
        start = ti(1); % ncread uses 1-based indices
        count = length(ti);
        Stokes_tmp = mean(ncread(dts_nc_in,'Stokes',[1 start],[nx count]),2);
        AntiStokes_tmp = mean(ncread(dts_nc_in,'AntiStokes',[1 start],[nx count]),2);
        tref_int_tmp = mean(ncread(dts_nc_in,'tref_int',start,count));
        tref_1_tmp = mean(ncread(dts_nc_in,'tref_1',start,count));
        tref_2_tmp = mean(ncread(dts_nc_in,'tref_2',start,count));
    end
    
    ncwrite(dts_nc_out,'Stokes',Stokes_tmp,[1 ii],[1 1]);
    ncwrite(dts_nc_out,'AntiStokes',AntiStokes_tmp,[1 ii],[1 1]);
    ncwrite(dts_nc_out,'tref_int',tref_int_tmp,ii,1);
    ncwrite(dts_nc_out,'tref_1',tref_1_tmp,ii,1);
    ncwrite(dts_nc_out,'tref_2',tref_2_tmp,ii,1);
end
toc
