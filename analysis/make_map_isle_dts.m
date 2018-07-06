%make_map_isle.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make_map_isle.m
%
%  %ISLE 2014 final locations of deployed moorings and DTS cables
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
run ../data_paths

position_file = [data_dir 'map_data/ISLE_asdeployed.mat'];
hires_file = [data_dir 'map_data/mvco_bathy.mat'];

warning off
load(position_file);
warning on

% Try to load high resolution bathy, download from OpenDAP server if
% necessary
hires = 1;
if hires & ~exist(hires_file,'file')
   try
       get_hires_bathy;
   catch
       hires = 0;
   end
end    
if hires
    load(hires_file)
end

config_name = 'asit';
lon_dts = ncread(dts_nc,'lon');
lat_dts = ncread(dts_nc,'lat');
d_dts = ncread(dts_nc,'distance');

% calibration points 
% (distances from calibrate_dts_linearsolve_gammathresh_towercal_cal6.m)
d3 = 143.0352;
d4 = 924.9522;
d5 = 3.9433e+03;
d6 = 4.5747e+03;

zi = find(~isnan(lon_dts));

lon3 = lon_dts(min(zi)); % First point is in coil before lat/lon points start
lat3 = lat_dts(min(zi));
lon4 = interp1(d_dts,lon_dts,d4);
lat4 = interp1(d_dts,lat_dts,d4);
lon5 = interp1(d_dts,lon_dts,d5);
lat5 = interp1(d_dts,lat_dts,d5);
lon6 = interp1(d_dts,lon_dts,d6);
lat6 = interp1(d_dts,lat_dts,d6);

c= load([data_dir 'map_data/coastline_MAB.mat']);
load([data_dir 'map_data/eastcoastbathy']);

[TLON,TLAT] = meshgrid(tlon,tlat);
tni = find(TLON < -76 & TLAT > 42);
tdata(tni) = NaN;

%%

figure(1); clf
landcolor=[245 222 179]./255;

subplot(2,1,1)
shelf_area = [-(74+30/60)  -(69+20/60)   40+30/60   42+10/60];
mvco_area = [-(70+38.5/60)  -(70+25/60)   41+15/60   41+21.2/60];

m_proj('Mercator','lon',shelf_area(1:2),'lat',shelf_area(3:4));

m_plot(c.lon_cst,c.lat_cst,'-','color','k','linewidth',1)

m_grid('tickdir','out');
hold on

%form color map that goes from white to off-blue
n=12;
bb=[.03 .46 .63];
bb2=[];
for i=1:3;
    bb2(:,i)=[bb(i):(1-bb(i))/n:1];
end
colormap(bb2)
m_contourf(tlon,tlat,tdata,[-1000:5:-5],'linestyle','none')
m_contour(tlon,tlat,tdata,[-1000:20:-20],'color',[.6 .6 .6])
caxis([-75 0])
m_plot([mvco_area(1), mvco_area(2), mvco_area(2), mvco_area(1), mvco_area(1)],....
       [mvco_area(3), mvco_area(3), mvco_area(4), mvco_area(4), mvco_area(3)],...
        '-','color','r','linewidth',0.5)
m_text(shelf_area(1)+5/60,shelf_area(4)-5/60,'a)','fontsize',8)    
m_plot(-74.0060,40.7128,'ko','markerfacecolor','b','markersize',4)
m_text(-74.0395,40.7128,{'New','York'},'color','k',...
        'horizontalalignment','right',...
        'verticalalignment','bottom','fontsize',8)
m_plot(-70.663184,41.526730,'ko','markerfacecolor','b','markersize',4)
m_text(-70.663184,41.536730,{'   Woods Hole'},'color','k',...
        'horizontalalignment','left',...
        'verticalalignment','top','fontsize',8)
    
subplot(2,1,2)
m_proj('Mercator','lon',mvco_area(1:2),'lat',mvco_area(3:4));

if hires
    m_contour(blon,blat,bh,[0,0],'k','linewidth',1);   
else
    m_plot(z(:,1),z(:,2),'k','linewidth',1.5)
end
    
m_grid('tickdir','out');
hold on
m_plot(ISLE_locs([4:5],1),ISLE_locs([4:5],2),'ks','markersize',5)
m_plot(ISLE_locs([3,6:end],1),ISLE_locs([3,6:end],2),'ks','markerfacecolor','k','markersize',5)
for ii=3:length(ISLE_loc_labels)
    m_text(ISLE_locs(ii,1)-0.7/60,ISLE_locs(ii,2)-0.2/60,ISLE_loc_labels(ii),'fontsize',8)
end
m_plot(BS_locs(:,1),BS_locs(:,2),'ks','markerfacecolor','k','markersize',5)

ii=1;    m_text(BS_locs(ii,1)-0.6/60,BS_locs(ii,2)-0.2/60,'E','fontsize',8)
ii=2;    m_text(BS_locs(ii,1)-0.6/60,BS_locs(ii,2),'I','fontsize',8)

%form color map that goes from white to off-blue
n=12;
bb=[.03 .46 .63];
bb2=[];
for i=1:3;
    bb2(:,i)=[bb(i):(1-bb(i))/n:1];
end
colormap(bb2)

if hires
    m_contourf(blon,blat,bh,[-50:5:-5],'color','none');
    m_contour(blon,blat,bh,[-80:5:-5],'color',[.6 .6 .6]);    
else 
    m_contourf(lonb,latb,bathb,[-50:5:-5],'color','none');
    m_contour(lonb,latb,bathb,[-80:5:-5],'color',[.6 .6 .6]);
end
m_text(mvco_area(1)-2/60,mvco_area(4)-0.1/60,'b)','fontsize',8)    
caxis([-75 0])

m_plot(lon3,lat3,'b.','markersize',16)
m_plot(lon4,lat4,'b.','markersize',16)
m_plot(lon5,lat5,'b.','markersize',16)
m_plot(lon6,lat6,'b.','markersize',16)

m_text(lon3+0.1/60,lat3+0.1/60,'R3/ASIT','color','b','fontsize',8)
m_text(lon4-0.4/60,lat4+0.3/60,'R4','color','b','fontsize',8)
m_text(lon5-0.1/60,lat5-0.4/60,'R5','color','b','fontsize',8)
m_text(lon6+0.2/60,lat6+0.2/60,'R6','color','b','fontsize',8)

m_plot(lon_dts,lat_dts,'r-','linewidth',2)
m_text(-70.4750,41.2917,{'Wasque','Shoals'},'color','k','fontsize',8,...
                'horizontalalignment','center')
            
scale_lat = 41+15.5/60;
scale_lon = -(70+38/60);
scale_len = 1/60;
m_quiver(scale_lon,scale_lat,0,scale_len,0,'k','linewidth',1)
m_quiver(scale_lon,scale_lat,scale_len,0,0,'k','linewidth',1)
m_text(scale_lon+scale_len+0.008,scale_lat,'x','fontsize',8,'horizontalalignment','left')
m_text(scale_lon,scale_lat+scale_len+0.003,'y','fontsize',8,'verticalalignment','bottom')

set(gcf,'renderer','zbuffer')
set(gcf, 'PaperSize', [4.5 5]); %Keep the paper size [width height] 
set(gcf, 'PaperPosition', [0 0.1 4.5 5]); %
print('-dpdf',['../figures/fig_ISLE_DTS_map.pdf'])