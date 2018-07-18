data_dir = '../../DTS-MVCO-data/';

microcat_dir = [data_dir 'ISLE_mat/microcats/'];
seagauge_dir = [data_dir 'ISLE_mat/seagauge/'];
wtpro_dir = [data_dir 'DTS_MVCO_cal/WaterTempPro/'];
adcp_dir = [data_dir 'ISLE_mat/ADCP/'];
adcp_nc_dir = [data_dir 'archive_Isle_mooring/'];
dts_dir = [data_dir 'DTS_MVCO_nc/'];
quadpod_dir = [data_dir 'ISLE_mat/quadpod/'];

dts_nc = [dts_dir 'DTSasit_proc.nc'];
cal_nc = [dts_dir 'DTSasit_cal.nc'];

microcat_mat = [microcat_dir 'microcats.mat'];
seagauge_mat = [seagauge_dir 'seagauge.mat'];

adcp_c_mat = [adcp_dir 'ISLE_C__burstmeans_06-Mar-2015.mat'];
adcp_g_mat = [adcp_dir 'ISLE_G__burstmeans_06-Mar-2015.mat'];
adcp_h_mat = [adcp_dir 'ISLE_H_corrected_11252015.mat'];
adcp_i_mat = [adcp_dir 'ISLE_I_burstmeans_06-Mar-2015.mat'];

cal_dir = [data_dir 'DTS_MVCO_cal/'];
wtcalfile = [cal_dir 'WaterTempPro/temppro_cals.mat'];

quadpod_e_mat = [quadpod_dir 'ISLE_E_quadpod_temp_ave_08172016.mat'];
