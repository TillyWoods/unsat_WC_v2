% Debugging runs - 14/10/2026

% Trying to identify why the Hnew = M\b equation in
% enthalpy_water_content_solve.m becomes singular at t = 40 (i.e. switch
% from freeze to melt phase) in all Kd=0 runs.

addpath(genpath('.'))   % Add all subfolders to path

% Load one of the runs that has the issue.
% load('Outputs\wc_periodic_switch_Kd_0_dZ_0-01_dt_0-005_14-Oct-2026_15_17_33.mat')
load('Outputs\14-10-2026\wc_periodic_switch_Kd_20_dZ_0-000625_dt_1-25e-05_15-Oct-2026_04_59_52.mat')
ss = ss_save;

% Save every timestep for debugging
dd.gap = 1;
dd.dt = dd.dt/2;

% Use the final successful timestep that is saved as the initial condition,
% remembering to also use the final time as the initial time (to get the
% periodic switching forcing correct)
ss0 = ss(end); %ss(end);
t0 = ss0.t;

tspan = [t0,t0+1];

% No need to calculate dd, pp etc., because these are loaded from the
% previous solution above.

% Timestepping
tic
ss = timestepping(ss0,pp,dd,tspan);
toc