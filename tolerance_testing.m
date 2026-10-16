function tolerance_testing(tolS1vals)

% Tolerance testing
% 16/10/2026

% For a fixed resolution (dZ), test the effect of changing the tolerances
% tolphi0 and tolS1.

tolphi0vals = [1e-3,1e-4,1e-5];
%tolS1vals = [1e-3,1e-4,1e-5];

% Lateral drainage values
Kdvals = [0,20];

% Time span and discretisation
tspan = [0,120]; %[0,300];
dZ = 0.005;


% Run each resolution with each drainage value
for kk = 1:length(Kdvals)
    Kd = Kdvals(kk);
    if Kd == 0
        dt = 0.4*0.5*1e-2;
    else
        dt = 0.5*0.25*1e-3;
    end
    for j = 1:length(tolS1vals)
        for i = 1:length(tolphi0vals)
            tic
            fprintf('Starting run: Kd = %i, tolphi0 = %f, tolS1 = %f', Kd, tolphi0vals(i), tolS1vals(j))
            periodic_switching_run(dt,dZ,Kd,tspan,tolphi0=tolphi0vals(i),tolS1=tolS1vals(j))
            fprintf('Done: Kd = %i, tolphi0 = %f, tolS1 = %f', Kd, tolphi0vals(i), tolS1vals(j))
            toc
        end
    end
end

end