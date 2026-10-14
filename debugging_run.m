% Debugging runs - 14/10/2026

% Trying to identify why the Hnew = M\b equation in
% enthalpy_water_content_solve.m becomes singular at t = 40 (i.e. switch
% from freeze to melt phase) in all Kd=0 runs.

% Load one of the runs that has the issue.
load('Outputs\wc_periodic_switch_Kd_0_dZ_0-01_dt_0-005_14-Oct-2026_15_17_33.mat')
ss = ss_save;

% Use the final successful timestep that is saved as the initial condition,
% remembering to also use the final time as the initial time (to get the
% periodic switching forcing correct)
ss0 = ss(end);
t0 = ss0.t;

tspan = [t0,t0+60];

% No need to calculate dd, pp etc., because these are loaded from the
% previous solution above.

% Timestepping
tic
ss = timestepping(ss0,pp,dd,tspan);
toc