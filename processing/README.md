# Data Processing Steps

* Convert raw DTS data from xml to NetCDF: run make_asit_dts_file.m
* Average DTS data in NetCDF file: timeavg_asit_dts_ncfile.m (intermediate)
* Average calibration data: timeavg_combine_asit_dts_ncfile (intermediate)
* Calibrate to create processed dataset (intermediate)
