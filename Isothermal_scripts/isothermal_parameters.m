function pp = isothermal_parameters(pp)
    %% Dimensionless parameters to use for the isothermal water transport case
    % Modifies the default dimensionless parameters pp.
    pp.Qsi = @(t) 0;
    pp.Q0 = @(t) 0;
    q0 = 0.1;
    pp.Rr = @(t) q0;
    pp.K = 100;
    pp.B = 1e4;
    pp.Tinf = 0;
    pp.Kd = 10;     % subsurface drainage
    pp.tolS1 = 1e-10;
end