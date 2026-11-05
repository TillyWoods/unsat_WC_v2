%% Script to make periodic switching snapshot plots
% Figures 6 and 7 in the manuscript

output_dir = ".";

% Load data
load("wc_periodic_switch_Kd_20_dZ_0-00125_dt_2-5e-05_23-Oct-2026_22_49_27.mat")
ss = ss_save;

% Freezing snapshots (fig 6)
plot_times_freeze = [22.2,22.5,22.8,23,24,25]-2;
plotting_snapshots_periodic_switching_func(ss,pp,plot_times_freeze,fluxes=1,output_dir=output_dir)

% Melting snapshots (fig 7)
plot_times_melt = [43,44,44.4,45,47,50.5]-2;
plotting_snapshots_periodic_switching_func(ss,pp,plot_times_melt,output_dir=output_dir)