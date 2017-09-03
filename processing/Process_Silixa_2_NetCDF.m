function Process_Silixa_2_NetCDF(dirname,Dstart,Dstop,outfile,append)
% Process_Silixa Version 1 (7/7/2014) 
% Process_Silixa: fuction to parse the *xml data files created by Silixa
%           DTS (Distributed Temperature Sensor) systems, specifically
%           tested for output from XT and Ultima systems.  The routine
%           detects external temperature sensors, as well as whether the
%           measurements were taken in a single-ended or double-ended
%           configuration, and saves the file appropriately.
%
% USAGE: Process_Silixa(dirname,Dstart,Dstop,outfile)
%        dirname = name and path of the directoyr containing XML files
%            (make sure to include the trailing '\' (windows) or '/' (mac))
%        Dstart = beginning (in meters)
%        Dstop = and end (m), or min and max distances on the DTS trace to 
%            be processed
%        outfile = output filename
%
% Example: Process_Silixa('ddf_files/',0,2016,'DTS_Trial.mat')
% Mac Example: Process_Silixa('/Users/myname/Desktop/DTS_data/full data set/deployment name/channel 1/2009/jul/',0,2076,'Deployment_Jul09.mat')
% PC Example: Process_Silixa('C:\Username\DTS Projects\full data set\deployment name\channel 3\2009\sep\',-50,1045,'Deployment_Sep09.mat')
%
% INPUTS: (variables) dirname,Dstart,Dstop,outfile
%         (files)   *.xml
%
% OUTPUT: (file) outfile.mat, containing the following (variables):
%       datetime = date and time of each measurement
%       tref_int = internal refrence temperature (degrees C)
%    	distance = measurement points along fiber (in meters)
%      	tempC = temperature of fiber (degrees C)
%    	Stokes = raw Stokes (dimensionless) amplitude along fiber
%      	AntiStokes = raw AntiStokes (dimensionless) amplitude along fiber
% OUTPUT: (file) Depending on the data files, the output datafile may also
%         include the following variables:
%    	StokesR = raw Stokes (dimensionless) amplitude from the reverse
%           trace
%     	AntiStokesR = raw AntiStokes (dimensionless) amplitude from the
%          	reverse trace
%      	tref_1 = external refrence temperature (1) (degrees C)
%    	tref_2 = external refrence temperature (2) (degrees C)
%       datetime_offset = datetime offset from DTS configuration
%
% Process_Silixa has been formated to match the layout, design, and output
% of Process_Sensornet written by Mark Hausner (2011).
%
% WRITTEN BY:   Scott Kobs, 7/16/2014
%
% -----------------------
% Process_Silixa_2_NetCDF(dirname,Dstart,Dstop,outfile,append)
% Adapted from Process_Silixa by Tom Connolly, 12/30/2014
% to write to NetCDF file instead of mat-file

% Additional input option:
% append - 1 if append to existing NetCDF file, 0 if create new (default)


if nargin == 4
    append = 0;
elseif ~isnumeric(append)
    append = 0;
end

close all; clc
q = cputime;

if exist(outfile)~= 0
    warning(['OUTFILE ' outfile ' already exists'])
end

files=dir(strcat([dirname '*.xml']));
nf = length(files);

display(strcat(['Number of Data Files: ' num2str(nf)]));

datetime(1:nf) = 0;
tref_int(1:nf) = 0;

toffset = 0;
Ends = 0;
error_flag = 0;

for f = 1:nf

    xDoc = xmlread([dirname files(f).name]);

    l = xDoc.getElementsByTagName('minDateTimeIndex');
    temp = char(l.item(0).getFirstChild.getData);
    datetime(f) = datenum(temp(1:23),'yyyy-mm-ddTHH:MM:SS.FFF'); 
    if length(temp) > 24
        if f == 1
            toffset = 1;
            datetime_offset(1:nf) = 0;
        end
        datetime_offset(f) = str2num(temp(25:end));
    end

    l = xDoc.getElementsByTagName('referenceTemperature');
    tref_int(f) = str2num(l.item(0).getFirstChild.getData);
    
    l = xDoc.getElementsByTagName('probe1Temperature');
    probe_1 = str2num(l.item(0).getFirstChild.getData);
    if probe_1 ~= 850;
        if f == 1
            tref_1(1:nf) = 0;
        end
        tref_1(f) = probe_1;
    end
        
    l = xDoc.getElementsByTagName('probe2Temperature');
    probe_2 = str2num(l.item(0).getFirstChild.getData);
    if probe_2 ~= 850;
        if f == 1
            tref_2(1:nf) = 0;
        end
        tref_2(f) = probe_2;
    end
    
    l = xDoc.getElementsByTagName('logData');
    fiber = str2num(l.item(0).getTextContent);
    
    z = fiber(:,1);
    a = find((z >= Dstart).*(z <= Dstop));
    d = z(a);
    
    if f==1
        distance=d;       
        tempC(1:length(d),1:nf)=0;
        Stokes(1:length(d),1:nf)=0;
        AntiStokes(1:length(d),1:nf)=0;
        if size(fiber,2)>4;
            Ends = 2;
            StokesR(1:length(d),1:nf)=0;
            AntiStokesR(1:length(d),1:nf)=0;
        end
    end
    
    if length(distance)~=length(d)
        display('Error: distance vector changes.  Process only one configuration at a time.');
        error_flag=1;
        break
    end
    
    Stokes(:,f) = fiber(a,2);
    AntiStokes(:,f) = fiber(a,3);
    if Ends == 2
        StokesR(:,f) = fiber(a,4);
        AntiStokesR(:,f) = fiber(a,5);
        tempC(:,f) = fiber(a,6);
    else
        tempC(:,f) = fiber(a,4);
    end
    
    if mod(f,100)==0
        display(strcat([num2str(f) ': ' datestr(datetime(f))]));
    end
    
    clear xDoc l temp probe_1 probe_2 fiber z a d 
    
end

if error_flag>0
    x=Stokes(:,1:f-1); clear Stokes; Stokes=x; clear x;
    x=AntiStokes(:,1:f-1); clear AntiStokes; AntiStokes=x; clear x;
    x=tempC(:,1:f-1); clear tempC; tempC=x; clear x;
    x=datetime(1:f-1); clear datetime; datetime=x; clear x;
    x=tref_int(1:f-1); clear tref_int; tref_int=x; clear x;
    if exist('tref_1') == 1
        x=tref_1(1:f-1); clear tref_1; tref_1=x; clear x;
    end
   	if exist('tref_2') == 1
        x=tref_2(1:f-1); clear tref_2; tref_2=x; clear x;
    end
end

if ~issorted(datetime)
    warning('datetime not sorted - contact CTEMPs')
    
end

% -------------------- Write or append NetCDF file ----------------- %

nx = length(distance);
nt = length(datetime);

disp(['writing to file: ' outfile])
if append == 0
    
    % create new file, overwriting existing file
    ncid = netcdf.create(outfile,'CLOBBER');
    
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
    
    Temp_varID = netcdf.defVar(ncid,'tempC','NC_DOUBLE',[x_dimID, t_dimID]);
    netcdf.putAtt(ncid,Temp_varID,'variable','Temperature along cable (preliminary)')
    netcdf.putAtt(ncid,Temp_varID,'units','degrees Celcius')    
    
    tint_varID = netcdf.defVar(ncid,'tref_int','NC_DOUBLE',[t_dimID]);
    netcdf.putAtt(ncid,tint_varID,'variable','Internal Instrument Temperature')
    netcdf.putAtt(ncid,tint_varID,'units','degrees Celcius')  
    
    if exist('tref_1') == 1
        t1_varID = netcdf.defVar(ncid,'tref_1','NC_DOUBLE',[t_dimID]);
        netcdf.putAtt(ncid,t1_varID,'variable','Reference Temperature #1')
        netcdf.putAtt(ncid,t1_varID,'units','degrees Celcius')            
    end
    if exist('tref_2') == 1
        t2_varID = netcdf.defVar(ncid,'tref_2','NC_DOUBLE',[t_dimID]);
        netcdf.putAtt(ncid,t2_varID,'variable','Reference Temperature #2')
        netcdf.putAtt(ncid,t2_varID,'units','degrees Celcius')    
    end
    if exist('datetime_offset')==1
        doff_varID = netcdf.defVar(ncid,'datetime_offset','NC_DOUBLE',[t_dimID]);
    end
    % end define mode and enter data mode
    netcdf.endDef(ncid);
        
    % write variables    
    netcdf.putVar(ncid,dist_varID,0,nx,distance)
    netcdf.putVar(ncid,date_varID,0,nt,datetime)
    netcdf.putVar(ncid,AntiStokes_varID,[0 0],[nx nt],AntiStokes)
    netcdf.putVar(ncid,Stokes_varID,[0 0],[nx nt],Stokes)
    netcdf.putVar(ncid,Temp_varID,[0 0],[nx nt],tempC)
    netcdf.putVar(ncid,tint_varID,0,nt,tref_int)
    if exist('tref_1') == 1
        netcdf.putVar(ncid,t1_varID,0,nt,tref_1)
    end    
    if exist('tref_2') == 1
        netcdf.putVar(ncid,t2_varID,0,nt,tref_2)
    end    
    if exist('datetime_offset')==1
        netcdf.putVar(ncid,doff_varID,0,nt,datetime_offset)
    end
    
else 
    ncid = netcdf.open(outfile,'WRITE');
    
    % find index to start writing time dimension (st)
    t_dimID = netcdf.inqDimID(ncid,'t');
    [~, st] = netcdf.inqDim(ncid,t_dimID);
    
    % write data
    date_varID = netcdf.inqVarID(ncid,'datetime');
    netcdf.putVar(ncid,date_varID,st,nt,datetime)
    
    AntiStokes_varID = netcdf.inqVarID(ncid,'AntiStokes');
    netcdf.putVar(ncid,AntiStokes_varID,[0 st],[nx nt],AntiStokes)
    
    Stokes_varID = netcdf.inqVarID(ncid,'Stokes');
    netcdf.putVar(ncid,Stokes_varID,[0 st],[nx nt],Stokes)
    
    Temp_varID = netcdf.inqVarID(ncid,'tempC');
    netcdf.putVar(ncid,Temp_varID,[0 st],[nx nt],tempC)

    tint_varID = netcdf.inqVarID(ncid,'tref_int');    
    netcdf.putVar(ncid,tint_varID,st,nt,tref_int)
    
    if exist('tref_1') == 1
        t1_varID = netcdf.inqVarID(ncid,'tref_1'); 
        netcdf.putVar(ncid,t1_varID,st,nt,tref_1)
    end    
    if exist('tref_2') == 1
        t2_varID = netcdf.inqVarID(ncid,'tref_2');
        netcdf.putVar(ncid,t2_varID,st,nt,tref_2)
    end    
    if exist('datetime_offset')==1
        doff_varID = netcdf.inqVarID(ncid,'datetime_offset');
        netcdf.putVar(ncid,doff_varID,st,nt,datetime_offset)
    end
end

netcdf.close(ncid);

end