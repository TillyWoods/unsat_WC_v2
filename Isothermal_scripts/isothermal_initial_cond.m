function ss0 = isothermal_initial_cond(phi0,dd,pp)
%% Isothermal transport initial conditions
% Input phi0 = "default" to use defalt value below. Otherwise, phi0 is a
% function handle, phi0 = @(Z) phi0.
if or(isstring(phi0),ischar(phi0))
    phi_func = @(Z) max(0,0.2*(exp(-Z/0.5)-exp(-Z(end)/0.5)));
else
    phi_func = phi0;
end
ss0.phi = phi_func(dd.Z);
ss0.S = 0*dd.Z.^0;
ss0.T = 0*dd.Z.^0;
ss0.msurf = 0;
ss0.mint = 0*dd.Z.^0;
ss0.hdot = 0;
ss0.roff = 0;
%ss0.Fminus = 0*dd.Zedge.^0;
%ss0.Fplus = 0*dd.Zedge.^0;
ss0.F = 0*dd.Zedge.^0;
ss0.Zm = NaN;
[ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp);
ss0.Zf = NaN;
ss0.q = 0*dd.Zedge.^0;
%ss0.a = pp.a;
ss0.h = 0;
end