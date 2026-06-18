%% Test run on tethys

dZ = 0.00125;           % resolution near the surface
dt = 0.25*0.25*1e-3;    % timestep size
Kd = 20;                % lateral drainage parameter
tspan = [0,300];        % time span

periodic_switching_run(dt,dZ,Kd,tspan)