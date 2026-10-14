%% Periodic switching resoltion testing

% dZ and corresponding dt values
dZvals = [1,1/2,1/4,1/8,1/16]*0.01;
dtvals20 = [1,0.5,0.25,0.1,0.05]*0.25*1e-3;
dtvals0 = [1,0.4,0.15,0.06,0.03]*0.5*1e-2;

% Lateral drainage values
Kdvals = [20,0];

% Time span
tspan = [0,300];

% Run each resolution wuth each drainage value
for i = 1:length(dZvals)
    for j = 1:length(Kdvals)
        if Kdvals(j)==0
            dt = dtvals0(i);
        else
            dt = dtvals20(i);
        end
        tic
        periodic_switching_run(dt,dZvals(i),Kdvals(j),tspan)
        toc
    end
end