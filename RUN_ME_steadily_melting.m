%% Steadily melting experiments

% Add all relevant folders to path
addpath(genpath('./..'))

%% Discretisation
dt = 0.25*1e-3;
dd = discretise_isothermal("dt",dt);

%% Parameters and forcings
% Dimensional parameters
dim = dimensional_parameters();

% Steadily melting forcing
Qsi = 200;
Q0 = 50;
dim.Qsi = @(t) Qsi;
dim.Q0 = @(t) Q0;

scale = scalings(dim);
pp = dimensionless_parameters(dim,scale);
pp.tolphi0 = 1e-4;

% Subsurface drainage?
pp.Kd = 3;      % values from 0 to 200 used in the manuscript

%% Initial condition
% Use steadily melting solution WITHOUT subsurface drainage as initial
% condition
ss0 = steadily_melting_initial_cond(pp.Qsi(1),pp.Q0(1),pp,dd,"drainage",0);
[ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp);

tspan = [0,50];

%% Timestepping
tic
ss = timestepping(ss0,pp,dd,tspan);
toc

%% Save data
ss_save = ss;
save("Outputs/wc_sm_Kd_"+num2str(pp.Kd,2)+"_"+strrep(strrep(string(datetime),' ','_'),':','_')+".mat","ss_save","pp","dd","scale","dim")