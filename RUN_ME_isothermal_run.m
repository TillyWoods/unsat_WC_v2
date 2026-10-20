%% Isothermal water transport run

% Add all relevant folders to path
addpath(genpath('./..'))

%% Discretisation
dt = 0.25*1e-2;
dd = discretise_isothermal("dt",dt);

%% Parameters and forcings
% Dimensional parameters
dim = dimensional_parameters();

% Dimensionless parameters
scale = scalings(dim);
pp = dimensionless_parameters(dim,scale);
pp = isothermal_parameters(pp);     % update to isothermal values

% Subsurface drainage?
pp.Kd = 3;     % values of 0, 2, 10 used in manuscript

%% Initial condition
phi00 = 0.2;
Z0 = 0.5;
phi0 = @(Z) phi00*exp(-Z/Z0);

% Create necessary structure of the initial condition
ss0 = isothermal_initial_cond(phi0,dd,pp);

% Time span
tspan = [0,5];

%% Timestepping
tic
ss = timestepping(ss0,pp,dd,tspan);
toc

%% Save data
ss_save = ss;
save("Ouputs/wc_isothermal_Kd_"+num2str(pp.Kd,2)+"_"+strrep(strrep(string(datetime),' ','_'),':','_')+".mat","ss_save","pp","dd","scale","dim")