%ISLE_mooring_locs.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ISLE_mooring_locs.m
%
%  %ISLE 2014 final locations of deployed moorings and DTS cables
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
load_data_paths;

warning off
position_file = '../data/ISLE_asdeployed.mat';
load(position_file);
warning on

config_name = 'asit';
cal_nc = [dts_dir 'DTS' config_name '_cal.nc'];
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

c= load([data_dir 'coastline/MAB.mat']);
load([data_dir 'bathy/eastcoastbathy']);
b = load([data_dir 'bathy/crm_capecod.mat']);

[b.X,b.Y] = meshgrid(b.x,b.y);

[TLON,TLAT] = meshgrid(tlon,tlat);
tni = find(TLON < -76 & TLAT > 42);
tdata(tni) = NaN;

%%

figure(1); clf
landcolor=[245 222 179]./255;

subplot(2,1,1)
shelf_area = [-(74+30/60)  -(69+20/60)   40+15/60   42+30/60];
mvco_area = [-(70+44/60)  -(70+25/60)   41+15/60   41+22/60];

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
m_text(shelf_area(1)+10/60,shelf_area(4)-10/60,'a','fontsize',12)    
m_plot(-74.0060,40.7128,'ko','markerfacecolor','k','markersize',4)
m_text(-74.0060,40.7128,{'New','York'},'color','k',...
        'horizontalalignment','right',...
        'verticalalignment','bottom','fontsize',8)
m_plot(-71.0589,42.3601,'ko','markerfacecolor','k','markersize',4)
m_text(-71.0589,42.3601,{'Boston '},'color','k',...
        'horizontalalignment','right','fontsize',8)

    
subplot(2,1,2)
m_proj('Mercator','lon',mvco_area(1:2),'lat',mvco_area(3:4));
m_plot(z(:,1),z(:,2),'k','linewidth',1.5)
m_grid('tickdir','out');
hold on
m_plot(ISLE_locs(:,1),ISLE_locs(:,2),'ks','markerfacecolor','k','markersize',5)
for ii=1:length(ISLE_loc_labels)
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
m_contourf(lonb,latb,bathb,[-50 -45 -40 -35 -30 -25 -20 -15 -10 -5],'color','none');
m_contour(lonb,latb,bathb,[-80 -40 -20],'color',[.6 .6 .6]);
m_text(mvco_area(1)+0.6/60,mvco_area(4)-0.6/60,'b','fontsize',12)    
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

%%

set(gcf,'renderer','zbuffer')
set(gcf, 'PaperSize', [4.5 5]); %Keep the paper size [width height] 
set(gcf, 'PaperPosition', [0 0.1 4.5 5]); %
print('-dpdf',['figures_paper/fig_ISLE_DTS_map.pdf'])