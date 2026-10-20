%% Periodic switching experiments - 30/05/2025
%
% A periodic version of the switching experiment in the unsaturated
% chapter of my thesis.
%
% Using the implicit, unsaturated code (with lateral drainage), but for
% runs both with and without lateral drainage.
%
% The forcing goes periodically between 3 sections of constant forcing
% values, detailed below.

%%
addpath(genpath('.'))   % Add all subfolders to path

Kd_vals = 20; %[20,0];
dZf_vals = 0.0025; %0.0025; %0.00125; %[0.005,0.0025,0.00125,0.0025,0.00125];
dtmultiplier = 0.25; %0.5; %0.25; %[0.5,0.5,0.5,0.25,0.25]; %[0.5*0.25*1e-3,0.5*0.25*1e-3,0.5*0.25*1e-3,0.25*0.25*1e-3,0.25*0.25*1e-3];
dZdtname = ["dZ_0-0025_dtmult_0-25"]; %["dZ_0-00125_dtmult_0.125"]; % ["dZ_0-005_dt_small","dZ_0-0025_dt_small","dZ_0-00125_dt_small","dZ_0-0025_dt_tiny","dZ_0-00125_dt_tiny"];

% Use no drainage solution after 10 periods as initial condition, if the
% logical below is set to 1 (added 14/08/2025). Otherwise, use steadily
% melting initial condition
%load("wc_periodic_switch_Kd_0_14-Aug-2025_15_59_12.mat")
% load("wc_periodic_switch_Kd_3_19-Aug-2025_13_04_20.mat")   % 3 drainage periods run after the 10 no-drainage periods above (19/08/2025)
% load("wc_periodic_switch_Kd_0_19-Aug-2025_11_33_55.mat")   % 3 no drainage periods run after the 10 no-drainage periods above (19/08/2025)
% ss_init = ss_save(end);
use_10period_init = 0;


for i = 1:length(Kd_vals)
    for j = 1:length(dZf_vals)

        %%% Saturated or unsaturated?
        Kd = Kd_vals(i);

        runname = "Kd_"+num2str(Kd)+"_"+dZdtname(j);

        % % Use smaller timesteps if water flow is allowed
        % if Kd == 0
        %     dt = 0.5*0.25*1e-3; % 1e-2;
        %     gap = 2*4*1e2; %10;
        % else
        %     dt = 0.5*0.25*1e-3;
        %     gap = 2*4*1e2;
        % end

        dt = 0.25*1e-3*dtmultiplier(j);
        gap = 400/dtmultiplier(j);
        
        %%% Discretisation
        Zend = 100;         % size of domain
        Zvfine = 1.5; %2; %1.5;       % size of very fine domain near the surface
        dZedgevF = dZf_vals(j); %0.0025;  % grid size in very fine part of domain
        Zfine = 4;          % bottom of the fine domain that extends below the very fine
        dZedgeF = 0.01;     % grid size is fine part of domain
        dZedgeC = 1e-1;     % grid size in coarse part of domain
        NT = 10;            % number of grid cells in which to linearly transition between fine and coarse mesh
        Zedge = nonuniform_Zedge_3_parts(Zend,Zvfine,dZedgevF,Zfine,dZedgeF,dZedgeC,NT);
        
        dd = discretise(Zedge);
        
        %dd = discretise_nonuniform(); %discretise_isothermal();
        dd.plot = 0;
        dd.gap = gap; %1e2; %10;
        dd.dt = dt;
        dd.steady_stop = 0;
        dd.steady_tol = 1e-4; %1e-4;
        
        %%% Parameters and forcings
        dim = dimensional_parameters();
        % Values used for steadily melting initial condition
        Qsi = 100; %200; %0;
        Q0 = 50; %-20; %20; %100; 
        
        % First section:
        Qsi1 = 100;
        Q01 = 50;
        gap1 = 20;      % duration of section
        
        % Second section:
        Qsi2 = 0;
        Q02 = -50;
        gap2 = 20;      % duration of section
        
        % Third section:
        Qsi3 = 0;
        Q03 = 100;
        gap3 = 20;      % duration of section
        
        period = gap1 + gap2 + gap3;
        
        dim.Qsi = @(t) Qsi1*(mod(t,period*dim.td)<gap1*dim.td)...
            +Qsi2*and(mod(t,period*dim.td)>=gap1*dim.td,mod(t,period*dim.td)<=(gap1+gap2)*dim.td)...
            +Qsi3*(mod(t,period*dim.td)>(gap1+gap2)*dim.td);
        dim.Q0 = @(t) Q01*(mod(t,period*dim.td)<gap1*dim.td)...
            +Q02*and(mod(t,period*dim.td)>=gap1*dim.td,mod(t,period*dim.td)<=(gap1+gap2)*dim.td)...
            +Q03*(mod(t,period*dim.td)>(gap1+gap2)*dim.td);
        
        scale = scalings(dim);
        pp = dimensionless_parameters(dim,scale);
        
        % Drainage or not?
        pp.Kd = Kd;
        
        % pp.B = 1e4;
        
        pp.tolphi0 = 1e-4;
        pp.tolS1 = 1e-3;
        
        %%% Initial condition
        if use_10period_init==1
            % Use end of 10 period no-drainage run as initial condition (added
            % 14/08/2025)
            ss0 = ss_init;
        else
            % Use steadily melting as initial condition
            ss0 = steadily_melting_initial_cond(Qsi/scale.F,Q0/scale.F,pp,dd);  % Dimensionless Qsi and Q0 in here
            % ss0.phi = 0*dd.Z.^0;
            % ss0.S = 1*dd.Z.^0; %0.5*dd.Z.^0;        % fully saturated initially
            [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp);
        end
        
        tspan = [0,5*period]; %[0,20*period]; %[0,10*period]; %[0,1]; %[0,100]; %[0,0.1];
        
        
        %% Timestepping
        tic
        ss = timestepping(ss0,pp,dd,tspan); %,runname);
        toc
        
        %% Save data
        save_time = 240; %inf;      % save the last save_time amount of time of data
        save_ind = ss(end).t-[ss.t]<save_time;
        % save_ind = 1:40:length(ss);
        ss_save = ss(save_ind);
        tsave = [ss_save.t];
        
        save(["wc_periodic_switch_Kd_"+num2str(pp.Kd,2)+"_"+strrep(strrep(string(datetime),' ','_'),':','_')+".mat"],"ss_save","pp","dd","scale","dim")

    end
end