% detide_g
% --------
% Use t_tide to detide velocities from ISLE G mooring

const={'M2';'M6';'M8';'N2';'K1';'O1';'S1';'S2'};
dt=.5;

ij=find(G.M.mtime>datenum(2014,7,10) &  G.M.mtime<datenum(2015,2,10)) ;
tt=G.M.mtime(ij);
uu=G.M.evm(:,ij); vv=G.M.nvm(:,ij); 
utide=nan.*uu; vtide=nan.*uu;

 
for ii=1:length(G.M.z)
    u=uu(ii,:); v=vv(ii,:);
    ik=find(isnan(u)==0);
    if length(ik)>500
        a=ik(1):ik(end);
        [name,freq,tidecon,tout]=t_tide(u(a)+sqrt(-1).*v(a),'interval',dt,'latitude',41,'start time',[tt(a(1))],'rayleigh',const);
        utide(ii,a)=real(tout);
        vtide(ii,a)=imag(tout);
    end
end

G.utide = nan(size(G.M.evm));
G.vtide = nan(size(G.M.nvm));
G.utide(:,ij) = utide;
G.vtide(:,ij) = vtide;

