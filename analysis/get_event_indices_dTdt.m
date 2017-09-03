function [eventi,event_daten] = get_event_indices_dTdt(T,daten,thresh)

% find all times when temperature difference is below a critical threshold
% (strong negative gradient, cooler towards east)
% define threshold as mean - sd*(std dev)
% default sd = 2

%tdiff = t1-t2;
%thresh = nanmean(tdiff)-sd*nanstd(tdiff);
%thresh = nanmean(tdiff)-2*nanstd(tdiff);

dt = (daten(2)-daten(1))*24;
dTdt = NaN*T;
dTdt(2:end-1) = (T(3:end)-T(1:end-2))/dt;

ti = find(dTdt < thresh);

% find indices of event start
starti = ti(find(diff([0; ti]) > 1));
endi = ti(find(diff([ti; length(dTdt)]) > 1));
if length(endi) == length(starti)-1
    endi = [endi;length(dTdt)];
end

% find indices of minimum (most negative) temperature difference
eventi = nan(size(starti));
for ii = 1:length(starti)
    alli = starti(ii):endi(ii);
    [~,mini] = min(dTdt(alli));
    eventi(ii) = alli(mini);
end
event_daten = daten(eventi);

% ensure that all dates are separated by a minimum time difference,
% to avoid double-counting events
mindt = 0.25; % 6 hours
sepi = diff(event_daten) > mindt;

eventi = eventi(sepi);
event_daten = event_daten(sepi);