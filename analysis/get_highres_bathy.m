% script get_highres_bathy.m
% --------------------------
% Download high-resolution bathymetry from NGDC 

disp('downloading MVCO bathy subset from NGDC OpenDAP server...')

lat_range = [41.2, 41.4];
lon_range = [-70.75, -70.4];

bathy_file_nc = 'https://www.ngdc.noaa.gov/thredds/dodsC/regional/nantucket_13_mhw_2008.nc';

lat_all = nc_varget(bathy_file_nc,'lat');
lon_all = nc_varget(bathy_file_nc,'lon');

ii = find(lat_all >= lat_range(1) & lat_all <= lat_range(2));
jj = find(lon_all >= lon_range(1) & lon_all <= lon_range(2));

bh = nc_varget(bathy_file_nc,'Band1',[ii(1),jj(1)],[length(ii),length(jj)]);

blat = lat_all(ii);
blon = lon_all(jj);

save ../data/mvco_bathy blon blat bh

disp('... bathy file saved to data/mvco_bathy.mat')