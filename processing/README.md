# Data Processing Steps

* Create averaged DTS data in NetCDF file: timeavg_asit_dts_ncfile.m 
* Average calibration data: timeavg_combine_asit_dts_ncfile 
* Calibrate to create processed, quality-controlled dataset: calibrate_dts_linearsolve_gammathresh_distavg_writenc.m

Note: These scripts process DTS data that has already been converted to NetCDF format. Software for processing raw XML can be found on the CTEMPS website at http://ctemps.org/data-processing