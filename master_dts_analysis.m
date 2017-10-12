% master_dts_analysis.m
% ---------------------
% Master script for analyzing DTS data and reproducing results

if ~exist('figures','dir')
    mkdir('figures')
end

cd analysis

make_map_isle_dts
evaluate_dts_calibration
timeseries_dts_dTdx_deployment
timeseries_dts_moorings_july
timeseries_dts_moorings_sept
event_m2_phase
event_composite_h
event_individual_phase_c

cd ..