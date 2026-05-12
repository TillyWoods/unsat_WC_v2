function pp = dimensionless_parameters(dim,scale)
    % 27th October 2022
    % Structure pp of dimensionless parameters and forcings, calculated
    % from the structures of dimensional parameters dim and scalings scale.
    % The forcings Qsi, Q0, Rs and Rr are functions of dimensionless time.
    
    pp.K = dim.k0*dim.rho*dim.g*scale.t/(dim.mu*scale.Z);       % dimensionless measure of permeability
    pp.khat = dim.k*scale.t/(dim.rho*dim.c*scale.Z^2);          % dimensionless thermal conductivity
    pp.Upsilon = scale.t*scale.F/(dim.rho*dim.c*scale.Z*scale.T); % dim.nu*scale.t/(dim.rho*dim.c*scale.Z);        % dimensionless heat transfer coefficient
    pp.B = dim.rho*dim.g*scale.Z*dim.dp/dim.gamma;              % Bond number (effect of gravity vs. cappillary pressure)
    pp.Stefan = dim.L/(dim.c*scale.T);                          % Stefan number (ration of latent heat to sensible heat)
    pp.a = dim.a;                                               % albedo
    pp.chi = dim.chi;                                           % proportion of absorbed shortwave radiation absorbed at the surface
    pp.kappa = dim.kappa*scale.Z;                               % extinction coefficient
    pp.s = dim.s;                                               % unmnamed optical property
    pp.alpha0 = dim.alpha0*scale.Z;
    pp.r0 = dim.r0*scale.Z;
     
    % pp.Kd = dim.Kd;                                             % lateral drainage scale/size
    pp.Ld = dim.Ld/scale.Z;                                      % length of lateral drainage catchment relative to vertical lengthscale
    pp.Kd = dim.Kd*scale.Z/dim.Ld;                               % dimensionless lateral drainage parameter (surface slope relative to horizontal and vertical lengthscales)

    pp.perm = @(phi) dim.perm(phi)/dim.k0;                       % permeability function k(phi)
    pp.krw = @(S) dim.krw(S);                                   % relative permeability function krw(S)
    
    pp.tol = dim.tol;         % set 1/phi = 1/max(phi,tol) to prevent issues when phi is small
    pp.tolphi0 = dim.tolphi0;       % (13/5/24) set phi = 0 when phi <= tolphi
    pp.tolS0 = dim.tolS0;           % (13/5/24) set S = 0 when S <= tolS0
    pp.tolS1 = dim.tolS1;           % (13/5/24) set S = 1 when S >= 1-tolS1

    % Microbes and nutrients (22/02/2023)
    pp.kappaF = dim.kPAR./(dim.alphaPAR.*scale.F);
    pp.kappaC = dim.kC/scale.Cw;
    pp.GA = dim.betaA.*scale.t;
    pp.GC = dim.betaC.*scale.t.*scale.Aw./scale.Cw;
    pp.deltaw = dim.dw.*scale.t;
    pp.deltai = dim.di.*scale.t;
    pp.DA = dim.DA.*scale.t./scale.Z.^2;
    pp.DC = dim.DC.*scale.t./scale.Z.^2;
    
    pp.Aiinf = dim.Aiinf/scale.Ai;
    pp.Ciinf = dim.Ciinf./scale.Ci;
    
    %% Forcings - most as functions of dimensionless time
    pp.Tinf = dim.Tinf/scale.T;                     % ice temperature at depth
    pp.Qsi = @(t) dim.Qsi(scale.t*t)/scale.F;       % shortwave radiation at surface
    pp.Q0 = @(t) dim.Q0(scale.t*t)/scale.F;         % combined longwave radiation and turbulent heat flux at surface
    pp.Rs = @(t) dim.Rs(scale.t*t)/scale.Rs;        % snowfall rate
    pp.Rr = @(t) dim.Rs(scale.t*t)/scale.Rr;        % rainfall rate
    pp.alpha = @(phi,S,t) dim.alpha(phi,S,t*scale.t)*scale.Z;               % absorption coefficient as a function of porosity and saturation
    pp.r = @(phi,S,t) dim.r(phi,S,t*scale.t)*scale.Z;                              % scattering coefficient as a function of porosity and saturation
    pp.nuhat = @(t)dim.nu(t*scale.t)*scale.T/scale.F;       % effective heat transfer coefficient
end