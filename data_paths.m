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

cal_dir = [data_dir 'DTS_MVCO_cal/'];
wtcalfile = [cal_dir 'WaterTempPro/temppro_cals.mat'];

quadpod_e_mat = [quadpod_dir 'ISLE_E_quadpod_temp_ave_08172016.mat'];
