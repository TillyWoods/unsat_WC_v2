% Restart test
% 23/10/2026

% To run some simulations for longer, because I think they might not have 
% reached a steady state.

addpath(genpath('.'))   % Add all subfolders to path

% Load the run to run further.
% load('.\..\Outputs\14-10-2026\wc_periodic_switch_Kd_20_dZ_0-00125_dt_2-5e-05_14-Oct-2026_23_57_24.mat')
load('.\..\Outputs\14-10-2026\wc_periodic_switch_Kd_20_dZ_0-01_dt_0-00025_14-Oct-2026_15_37_15.mat')
ss = ss_save;


% Use the final timestep that is saved as the initial condition,
% remembering to also use the final time as the initial time (to get the
% periodic switching forcing correct)
ss0 = ss(end); %ss(end);
t0 = ss0.t;

tspan = [t0,t0+180];

% No need to calculate dd, pp etc., because these are loaded from the
% previous solution above.

% Timestepping
tic
ss = timestepping(ss0,pp,dd,tspan);
toc