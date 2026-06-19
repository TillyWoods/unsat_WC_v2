%% Test run on tethys

dZ = 0.0025; %0.00125;           % resolution near the surface
dt = 0.2*0.5*1e-2; %0.1*0.25*1e-3;    % timestep size
Kd = 0;                % lateral drainage parameter
tspan = [0,300];        % time span

tic
periodic_switching_run(dt,dZ,Kd,tspan)
toc