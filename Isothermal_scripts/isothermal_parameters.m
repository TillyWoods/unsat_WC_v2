function pp = isothermal_parameters(pp)
% Conditions to use for the isothermal water transport case
% 21st Feb 2024
% Modifies the default dimensionless parameters pp.
pp.Qsi = @(t) 0;
pp.Q0 = @(t) 0;
q0 = 0.1;
pp.Rr = @(t) q0;
pp.K = 100;
pp.B = 100;
pp.Tinf = 0;
end