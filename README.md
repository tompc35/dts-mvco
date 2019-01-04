# DTS-MVCO

MATLAB code for processing, calibrating, evaluating and analyzing data from distributed temperature sensing (DTS) system deployed during the Inner Shelf Lateral Exchange (ISLE) project south of Martha's Vineyard in 2014.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

##### Reference

* Connolly, T. P. and A. R. Kirincich (2019) High-resolution observations of subsurface fronts and alongshore bottom temperature variability over the inner shelf, Journal of Geophysical Research. doi:[10.1029/2018jc014454](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018JC014454)

### Running the code
* *Optional* - Edit [data_paths.m](data_paths.m) so that `data_dir` variable points to the right location on your machine. This step is not necessary if the data and code directories are located in the same overarching directory:

```
data_dir/
  DTS-MVCO/
    analysis/
    processing/
    data_paths.m
    master_dts_analysis.m
    master_dts_processing.m
    README.md
  DTS-MVCO-data/
    archive_Isle_mooring/
    DTS_MVCO_cal/
    DTS_MVCO_nc/
```

* *Optional* - Run [master_dts_processing.m](master_dts_processing.m) to generate processed and calibrated DTS data in NetCDF format.
* Run [master_dts_analysis.m](master_dts_analysis.m) to generate figures.

### Data sources

DTS and associated calibration data are available on Zenodo:
* Connolly, T. P., and A. R. Kirincich (2018), Distributed temperature sensing and associated data - Martha’s Vineyard Coastal Observatory 2014, Zenodo. DOI:10.5281/zenodo.1136113 https://zenodo.org/record/1136113

ISLE mooring data are available on the WHOI Open Access Data Server:
* Kirincich, A. R., and S. J. Lentz (2017b), Inner shelf lateral exchange, Woods Hole Open Access Server (WHOAS), DOI:10.1575/1912/8740, https://hdl.handle.net/1912/8740


### Dependencies

* SNCTOOLS
-r4053
http://mexcdf.sourceforge.net

* t_tide
Version 1.02
https://www.eoas.ubc.ca/~rich/#T_Tide

* m_map
Version 1.4h
https://www.eoas.ubc.ca/~rich/map.html

* CSIRO seawater toolbox
Version 3.3.1
http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm

* RPSstuff
https://woodshole.er.usgs.gov/operations/sea-mat/RPSstuff-html/index.html

* movingstd
https://www.mathworks.com/matlabcentral/fileexchange/9428-movingstd---movingstd2

Code has been run using MATLAB R2015a (Mac) using versions of software specified above.
