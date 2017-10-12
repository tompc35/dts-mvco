% detide_c
% --------
% Use t_tide to detide velocities from ISLE C mooring
 
const={'M2';'M6';'M8';'N2';'K1';'O1';'S1';'S2'};
dt=.5;

 
%% 

ij=find(C.M.mtime>datenum(2014,7,10) &  C.M.mtime<datenum(2015,2,10)) ;
tt=C.M.mtime(ij);
uu=C.M.evm(:,ij); vv=C.M.nvm(:,ij); 
 
%detide all this data and save
utide=nan.*uu; vtide=nan.*uu;

 
for ii=1:length(C.M.z)
    u=uu(ii,:); v=vv(ii,:);
    ik=find(isnan(u)==0);
    if length(ik)>500
        a=ik(1):ik(end);
        [name,freq,tidecon,tout]=t_tide(u(a)+sqrt(-1).*v(a),'interval',dt,'latitude',41,'start time',[tt(a(1))],'rayleigh',const);
        utide(ii,a)=real(tout);
        vtide(ii,a)=imag(tout);
    end
end

C.utide = nan(size(C.M.evm));
C.vtide = nan(size(C.M.nvm));
C.utide(:,ij) = utide;
C.vtide(:,ij) = vtide;

