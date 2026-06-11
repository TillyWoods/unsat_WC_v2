function scale = scalings(dim)
    scale=struct;
    scale.Z = 1;                                        % length scale (m)
    scale.t = dim.td;                                   % timescale (s) eg. 1 day or 1 year in seconds
    scale.F = 200;                                      % radiation scale (W m^-2)
    scale.T = scale.t*scale.F/(dim.rho*dim.c*scale.Z); % scale.F/dim.nu;                           % temperature scale (K)
    scale.msurf = scale.Z/scale.t;                      % melt rate scale (m s^-1)
    scale.hdot = scale.msurf;                           % rate of change of surface height scale (m s^-1)
    scale.roff = scale.msurf;                           % surface runoff scale (m s^-1)
    scale.Rs = scale.msurf;                             % snowfall rate (m s^-1)
    scale.Rr = scale.msurf;                             % rainfall rate (m s^-1)
    scale.H = dim.rho*dim.c*scale.T;                    % enthalpy scale (J m^-3 = N m^-2)
    scale.q = scale.Z/scale.t;                          % water flux scale (m s^-1)
    scale.pw = dim.rho*dim.g*scale.Z;                   % pressure scale (Pa)
    scale.mint = dim.rho/scale.t;                       % internal melting scale (kg m^-3 s^-1)
end