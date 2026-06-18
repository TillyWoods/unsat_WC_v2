%% Periodic switching resoltion testing

% dZ and corresponding dt values
dZvals = [1,1/2,1/4,1/8,1/8,1/16]*0.01;
dtvals = [1,0.5,0.25,0.1,0.125,0.05]*0.25*1e-3;

% Lateral drainage values
Kdvals = [20,0];

% Time span
tspan = [0,300];

% Run each resolution wuth each drainage value
for i = 1:length(dZvals)
    for j = length(Kdvals)
        tic
        periodic_switching_run(dtvals(i),dZvals(j),Kdvals(j),tspan)
        toc
    end
end