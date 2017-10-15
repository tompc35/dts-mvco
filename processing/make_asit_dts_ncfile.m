% make_asit_dts_ncfile

clear all
run ../data_paths

nc_out_dir = dts_dir;
xml_base_dir = [data_dir 'DTS_MVCO_xml/XT14008/temperature/'];

all_deployments = {'asit_01','asit_02c'};
channel = 1;

for ii = 1 :length(all_deployments)
    deployment_name = char(all_deployments(ii))
    if strcmp(deployment_name,'asit_01')
        subfolders = {'July'};
    elseif strcmp(deployment_name,'asit_02c')
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
    else
        error('invalid deployment name')
    end

    tic
    diary('asit_dts_ncfile.log')
    for kk = 1:length(subfolders)
        disp(['Processing ' deployment_name ' - ' char(subfolders(kk))])

        xml_dir = [xml_base_dir ...
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
end
diary off