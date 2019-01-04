% master_dts_analysis.m
% ---------------------
% Master script for analyzing DTS data and reproducing results

if ~exist('figures','dir')
    mkdir('figures')
end

cd analysis

make_map_isle_dts             % Figure 1
timeseries_stats_dts          % Figures 2 and 6
evaluate_dts_calibration      % Figure 3
timeseries_dts_moorings_july  % Figure 4
timeseries_dts_moorings_sept  % Figure 5
event_m2_phase                % Figure 7
event_composite_h             % Figure 8
event_individual_phase_c      % Figures 9,10,11

cd ..