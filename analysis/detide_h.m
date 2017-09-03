% modified from script emailed by Anthony
% 12.10.2015

%clear all; close all;

load_dts_isle_data;

%cd ~/data/ISLE/
%load rawdata/ADCP/H/asitisleburstmeans_02_42_V01.mat

 
const={'M2';'M6';'M8';'N2';'K1';'O1';'S1';'S2'};

 
%load ISLE_H_corrected_11252015.mat
dt=.5;

 
%% 

 
%Just pull out tide for  august and sept
ij=find(H.ttime>datenum(2014,8,1) &  H.ttime<datenum(2015,9,10)) ;
tt=H.ttime(ij);
hu=H.uu(:,ij); hv=H.vv(:,ij); 
% cu=cuu(:,i); cv=cvv(:,i); 

 
%detide all this data and save
hutide=nan.*hu; hvtide=nan.*hu;

 
for ii=1:length(H.zz)
    u=hu(ii,:); v=hv(ii,:);
    ik=find(isnan(u)==0);
    m2_pha(ii) = NaN;
    m2_amp(ii) = NaN;
    if isempty(ik)==0
        a=ik(1):ik(end);
        [name,freq,tidecon,tout]=t_tide(u(a)+sqrt(-1).*v(a),'interval',dt,'latitude',41,'start time',[tt(a(1))],'rayleigh',const);
        hutide(ii,a)=real(tout);
        hvtide(ii,a)=imag(tout);
        m2_pha(ii) = tidecon(5,end-1);
        m2_amp(ii) = tidecon(5,1);
    end
end

H.utide = nan(size(H.utide));
H.vtide = nan(size(H.vtide));
H.utide(:,ij) = hutide;
H.vtide(:,ij) = hvtide;

