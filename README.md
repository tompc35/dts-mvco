# DTS-MVCO

MATLAB code for processing, calibrating, evaluating and analyzing data from
distributed temperature sensing (DTS) system deployed during the Inner Shelf
Lateral Exchange (ISLE) project south of Martha's Vineyard in 2014.

## Running the code
* Edit [data_paths.m](data_paths.m) so that `data_dir` variable points to the right location on your machine. By default, the data and code directories are assumed to be in the same overarching directory.
* *Optional* - Run [master_dts_processing.m](master_dts_processing.m) to generate processed and calibrated DTS data in NetCDF format.
* Run [master_dts_analysis.m](master_dts_analysis.m) to generate figures.

## Dependencies

* SNCTOOLS
-r4053
http://mexcdf.sourceforge.net

* t_tide
Version 1.02
https://www.eoas.ubc.ca/~rich/#T_Tide

* m_map
Version 1.4h
https://www.eoas.ubc.ca/~rich/map.html

* Code has been run using MATLAB R2015a on a MacBook Pro with 2.5 GHz Core i7
processor.
