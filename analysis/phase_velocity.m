function [cp,phixy,ia,ib,ic] = phase_velocity(lona,lata,Ta,lonb,latb,Tb,lonc,latc,Tc,dt)

% Compute phase velocity following Lee (1961), based on time of minimum
% dT/dt (i.e., most rapid cooling).
%
% Computes phase velocity of one event based on longitude,latitude and 
% temperature at three sites. The event arrival time at each site corresponds
% to the minimum in dT/dt.
%
% INPUTS:
% lon* - longitude of site * (a,b,c)
% lat* - latitude of site * (a,b,c)
% T* - temperature time series at site * (a,b,c)
% dt - time interval in temperature time series, in seconds
%
% OUTPUTS:
% cp - phase velocity magnitude (m/s)
% phixy - propagation direction (positive CCW from east)
% ia,ib,ic - indices of arrival times

% geometry of array
xa = 0; ya = 0;
[xc,yc] = latlon2xy(latc,lonc,lata,lona);
[xb,yb] = latlon2xy(latb,lonb,lata,lona);

% angle of deviation of a-b line from x (eastward) axis, degrees
theta = atan2(yb,xb)*180/pi;

% length of triangle segments
length_ca = sqrt((xc-xa)^2 + (yc-ya)^2);
length_ab = sqrt((xa-xb)^2 + (ya-yb)^2);
length_bc = sqrt((xb-xc)^2 + (yb-yc)^2);

% angles of triangle (law of cosines)
ang_ca = acos((-length_ca^2 + length_ab^2 + length_bc^2) ...
    /(2*length_ab*length_bc))*180/pi;
ang_ae = acos((-length_ab^2 + length_bc^2 + length_ca^2) ...
    /(2*length_bc*length_ca))*180/pi;
ang_bc = acos((-length_bc^2 + length_ca^2 + length_ab^2) ...
    /(2*length_ca*length_ab))*180/pi;

dTdta = nan(size(Ta));
dTdta(2:end-1) = 0.5*(Ta(3:end)-Ta(1:end-2))/dt;

dTdtb = nan(size(Tb));
dTdtb(2:end-1) = 0.5*(Tb(3:end)-Tb(1:end-2))/dt;

dTdtc = nan(size(Tc));
dTdtc(2:end-1) = 0.5*(Tc(3:end)-Tc(1:end-2))/dt;

[~,ic] = min(dTdtc);
[~,ia] = min(dTdta);
[~,ib] = min(dTdtb);

cp_ca = 1000*length_ca/((ic-ia)*dt);
cp_ab = -1000*length_ab/((ia-ib)*dt);
cp_bc = 1000*length_bc/((ib-ic)*dt);

phi = atan((cp_ab/cp_ca)*cscd(ang_bc) - cotd(ang_bc));
phi2 = atan(cotd(ang_ca) - (cp_ab/cp_bc)*cscd(ang_ca));

% magnitude of phase velocity
cp1 = cp_ab*cos(phi);
cp2 = cp_ca*cos(ang_bc*pi/180 - phi);
cp3 = cp_bc*cos(phi+ang_ca*pi/180);

if isfinite(cp1)
    cp = cp1;
elseif isfinite(cp2)
    cp = cp2;
elseif isfinite(cp3)
    cp = cp3;
else
    cp = NaN;
end

if cp < 0
    cp = -cp;
    phi = phi - pi;
end

% angle of propagation relative to x-axis
phixy = theta*pi/180-phi;

end