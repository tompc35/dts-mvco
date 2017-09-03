% modified from script emailed by Anthony
% 12.10.2015

%clear all; close all;

%load_dts_isle_data;

%cd ~/data/ISLE/
%load rawdata/ADCP/H/asitisleburstmeans_02_42_V01.mat

 
const={'M2';'M6';'M8';'N2';'K1';'O1';'S1';'S2'};

 
%load ISLE_H_corrected_11252015.mat
dt=.5;

 
%% 

 
%Just pull out tide for  july, august and sept
ij=find(I.M.mtime>datenum(2014,7,10) &  I.M.mtime<datenum(2015,2,10)) ;
tt=I.M.mtime(ij);
uu=I.M.evm(:,ij); vv=I.M.nvm(:,ij); 
% cu=cuu(:,i); cv=cvv(:,i); 

 
%detide all this data and save
utide=nan.*uu; vtide=nan.*uu;

 
for ii=1:length(I.M.z)
    u=uu(ii,:); v=vv(ii,:);
    ik=find(isnan(u)==0);
    if length(ik)>500
        a=ik(1):ik(end);
        [name,freq,tidecon,tout]=t_tide(u(a)+sqrt(-1).*v(a),'interval',dt,'latitude',41,'start time',[tt(a(1))],'rayleigh',const);
        utide(ii,a)=real(tout);
        vtide(ii,a)=imag(tout);
    end
end

I.utide = nan(size(I.M.evm));
I.vtide = nan(size(I.M.nvm));
I.utide(:,ij) = utide;
I.vtide(:,ij) = vtide;

