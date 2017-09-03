% make_asit_dts_matfile

clear all
addpath ../../tools/DTS/ -end

nc_out_dir = '/media/tompc/data/DTS_nc/';
xml_base_dir = '/media/tompc/data/DTS_data/XT14008/temperature/';

% deployment_name = 'asit_01';
% channel = 1;
% subfolders = {'July'};

deployment_name = 'asit_02c';
channel = 1;
subfolders = {  'DTS1_07b'
                'DTS1_07c'
                'DTS1_07d'
                'DTS1_08a'
                'DTS1_08b'
                'DTS1_08c'
                'DTS1_08d'
                'DTS1_09a'
                'DTS1_09b'
                'DTS1_09c'
                'DTS1_10a'
                'DTS1_10b'
                'DTS1_10c'
                'DTS1_10d'};

tic
diary('asit_dts_ncfile.log')
for kk = 1:length(subfolders)
    disp(['Processing ' deployment_name ' - ' char(subfolders(kk))])
    
    xml_dir = ['/media/tompc/data/DTS_data/XT14008/temperature/'...
                deployment_name '/channel' num2str(channel) '/'...
                char(subfolders(kk)) '/'];

    nc_out = [nc_out_dir 'DTS' deployment_name '_channel' num2str(channel) '.nc'];
    if kk == 1
        append = 0;
    else
        append = 1;
    end
    Process_Silixa_2_NetCDF(xml_dir,0,4873,nc_out,append)
    toc
end
diary off