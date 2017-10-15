% master_dts_processing.m
% -----------------------
% Master script for processing DTS data


cd processing

make_asit_dts_ncfile
timeavg_combine_asit_dts_ncfile
timeavg_asit_calfile
calibrate_dts_linearsolve_gammathresh_distavg_writenc

cd ..