%% Periodic switching experiments

%%
addpath(genpath('.'))   % Add all subfolders to path
        
%% Discretisation
dt = 0.5*0.25*1e-3;     % timestep size
dZedgevF = 0.005;        % 'very fine' grid resolution near surface
dZvfine = 2;            % thickness of 'very fine' zone
dd = discretise_nonuniform("dt",dt,"dZedgevF",dZedgevF,"Zvfine",dZvfine);

%% Parameters and forcings
% Dimensional parameters
dim = dimensional_parameters();

% Periodic switching forcings
[dim.Qsi,dim.Q0] = periodic_switching_forcing(dim);

scale = scalings(dim);
pp = dimensionless_parameters(dim,scale);

% Drainage or not?
pp.Kd = 20;

% Solid and saturated tolerances
pp.tolphi0 = 1e-4;
pp.tolS1 = 1e-3;

%% Initial condition
% Use steadily melting as initial condition
ss0 = steadily_melting_initial_cond(pp.Qsi(1),pp.Q0(1),pp,dd);  % Dimensionless Qsi and Q0 in here
[ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp);

tspan = [0,5*period];


%% Timestepping
tic
ss = timestepping(ss0,pp,dd,tspan);
toc

%% Save data
save_time = 240; %inf;      % save the last save_time amount of time of data
save_ind = ss(end).t-[ss.t]<save_time;
ss_save = ss(save_ind);
tsave = [ss_save.t];

save(["Outputs/wc_periodic_switch_Kd_"+num2str(pp.Kd,2)+"_"+strrep(strrep(string(datetime),' ','_'),':','_')+".mat"],"ss_save","pp","dd","scale","dim")