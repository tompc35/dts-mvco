data_dir = '/Users/tomconnolly/work/Data/';

microcat_dir = [data_dir 'ISLE/microcats/'];
seagauge_dir = [data_dir 'ISLE/seagauge/'];
wtpro_dir = [data_dir 'DTS_cal/WaterTempPro/'];
adcp_dir = [data_dir 'ISLE/ADCP/'];
dts_dir = [data_dir 'DTS_nc/'];
isle_node_dir = [data_dir 'ISLE/node/'];
quadpod_dir = [data_dir 'ISLE/quadpod/'];

dts_nc = [dts_dir '/proc/cal_6/DTSasit_proc.nc'];
cal_nc = [dts_dir 'DTSasit_cal.nc'];

microcat_mat = [microcat_dir 'microcats.mat'];
seagauge_mat = [seagauge_dir 'seagauge.mat'];

adcp_c_mat = [adcp_dir 'ISLE_C__burstmeans_06-Mar-2015.mat'];
adcp_g_mat = [adcp_dir 'ISLE_G__burstmeans_06-Mar-2015.mat'];
adcp_h_mat = [adcp_dir 'ISLE_H_corrected_11252015.mat'];
adcp_i_mat = [adcp_dir 'ISLE_I_burstmeans_06-Mar-2015.mat'];

wtcalfile = [data_dir 'DTS_cal/WaterTempPro/temppro_cals.mat'];

node_mat = [isle_node_dir 'node'];

quadpod_e_mat = [quadpod_dir 'ISLE_E_quadpod_temp_ave_08172016.mat'];
