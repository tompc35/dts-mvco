% detide_h
% --------
% Use t_tide to detide velocities from ISLE H mooring

 
const={'M2';'M6';'M8';'N2';'K1';'O1';'S1';'S2'};
dt=.5;

 
%% 

ij=find(H.ttime>datenum(2014,8,1) &  H.ttime<datenum(2015,2,10)) ;
tt=H.ttime(ij);
hu=H.uu(:,ij); hv=H.vv(:,ij); 
 
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

