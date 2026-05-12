function [S,q,qd,d] = sm_drainage_ana_sat(Z,pp,Qsi,Q0)
if pp.Kd==0
    S = ones(size(Z));
    q = zeros(size(Z));
    qd = zeros(size(Z));
    d = 0;
else
% Steadily melting solutions
Tinf = pp.Tinf;
msurf = (pp.chi*(1-pp.a)*Qsi+Q0)/pp.Stefan;
hdot = -(Qsi*(1-pp.a)+Q0)./(pp.Stefan-Tinf);
phi0 = -(1-pp.chi)*(1-pp.a)*Qsi/(hdot*pp.Stefan);
Zm = -(1/pp.kappa)*log(hdot*Tinf/((1-pp.chi)*(1-pp.a)*Qsi));

phi_func = @(Z) phi0*(exp(-pp.kappa*Z)-exp(-pp.kappa*Zm)).*(Z<Zm);
dphidZ_func = @(Z) -pp.kappa*phi0*exp(-pp.kappa*Z);

% Solve for S in the saturated region
dSdZ = @(Z,S) ((-3*pp.K*phi_func(Z).^2.*S.^3+hdot*(1-S))*dphidZ_func(Z)...
    - pp.K*pp.Kd*phi_func(Z).^3*S.^3)./(3*pp.K*S.^2.*phi_func(Z).^3 + hdot*phi_func(Z));
Zstop = Zm - 0.01;   % stop before Zm where dSdZ blows up
S0_guess = 0.55; %0.2;
phisurf = phi_func(0);
surface_func = @(S0) pp.K*phisurf^3*S0.^3 + hdot*phisurf*S0 + hdot*(1-phisurf);
S0 = fsolve(surface_func,S0_guess);
S0(or(S0<0,S0>1)) = 1;      % set S0 to 1 if a non-physical value is given (02/09/2024)
[Zout,Sout] = ode15s(dSdZ,[0,Zstop],S0);

% Interpolate unsaturated S onto Z grid
S = interp1(Zout,Sout,Z);

phi = phi_func(Z);
W = 1-phi+phi.*S;       % W in unsat region

% Calculate unsaturated and saturated fluxes
qunsat = pp.K*phi.^3.*S.^3;

qsat = -pp.Kd*pp.K*phi0^3*(-(1/(3*pp.kappa))*(exp(-3*pp.kappa*Z)-exp(-3*pp.kappa*Zm))...
    +(3/2/pp.kappa)*(exp(-2*pp.kappa*Z)*exp(-pp.kappa*Zm)-exp(-3*pp.kappa*Zm))...
    -(3/pp.kappa)*(exp(-pp.kappa*Z)*exp(-2*pp.kappa*Zm)-exp(-3*pp.kappa*Zm))-exp(-3*pp.kappa*Zm)*(Z-Zm));

% The water table Zd is where unsat and sat q+hdot*W is continuous
id = find((qsat+hdot)<=(qunsat+hdot*W),1,'first');
Zd = Z(id);
Zd(isempty(Zd)) = 0;
d = Zm-Zd;      % thickness of saturated region

% Overall values of S and q (unsat above Zd, sat below)
Sall = S;
Sall(Z>Zd)=1;
qall = qunsat;
qall(Z>Zd) = qsat(Z>Zd);

S = Sall;
q = qall;

% Lateral drainage
qd = pp.Kd*pp.K*phi.^3.*S.^3;
end
end