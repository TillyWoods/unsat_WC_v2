% 27th June 2024

% Trying to explain why Zt tends to a constant in isothermal with drainage.

q0 = 0.1;
Z0 = 0.5;
phi0 = 0.2;
K = 100;
Kd = 3;
% Ztlim_func = @(Zt) q0 - (q0/Kd)*(1-exp(-Kd*Zt)) - (Z0/3)^2*Kd*K*phi0^3*exp(-3*Zt/Z0);
Ztlim_func = @(Zt) q0*exp(-Kd*Zt) - (Z0/3)*Kd*K*phi0^3*exp(-3*Zt/Z0);

Ztlim = fsolve(Ztlim_func,0.25)

Ztvals = linspace(0,1);
qout = Ztlim_func(Ztvals);
figure(1)
plot(Ztvals,qout)

% Ztlim = Ztvals(qout==max(*qout));