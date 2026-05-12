function dim = dimensional_parameters()
    dim = struct;
    dim.k = 2.1;        % kgms^-3K^-1 (thermal conductivity) [CP2010,MH2017 2.1]
    dim.Cstar = 0.002;  % dimless (transfer coefficient) [CP2010 0.002]
    dim.P = 80e3;       % Pa (atmospheric pressure) [CP2010 80e3]
    dim.u = 5;          % ms^-1 (wind speed) [CP2010 5]
    dim.epsilon = 0.97;    % dimless (surface emissivity) [MH2017 0.97, Liston1999 0.98, CP2010 0.94-0.99]
    dim.sigma = 5.7e-8;    % Wm^-2K^-4 (Stefan-Boltzmann const.)
    dim.Tm = 273.15;    % K (melting temp) [CP2010 273.15]
    dim.phi0 = 1;       % dimless (reference porosity)
    dim.L = 334000;     % m^2s^-2 (latent heat of melting) [MH2017 334000]
    dim.c = 2097;       % Jkg^-1K^-1 (specific heat capacity) [CP2010 2097 (of ice at 0oC), MH2017 2050]
    dim.a = 0.6;        % dimless (albedo) [MH2017 0.6]
    dim.rho = 910;      % kgm^-3 (ice density) [CP2010 910, MH2010 917]
    dim.chi = 0.36;     % dimless (proportion of absorbed SW radiation abosrbed at surface)
    dim.kappa = 1.5;    % m^-1 (extinction coefficient) [TF2005 1.5]
    dim.dp = 1e-4;      % m (typical grain size) [MH2017 10^-4 m]
    dim.gamma = 0.07;   % N m^-1 (surface tension) [MH2017 0.07 N m^-1]
    dim.g = 9.81;       % m s^-2 (acceleration due to gravity)
    dim.mu = 1e-3;      % Pa s (viscosity of water) [MH2017 10^-3 Pa s]
    dim.ty = 365*24*60*60;    % s (number of seconds in a year)
    dim.td = 24*60*60;        % s (number of seconds in a day)
    dim.th = 60*60;           % s (number of seconds in a hour)
    
    dim.tol = 1e-6;         % set 1/phi = 1/max(phi,tol) to prevent issues when phi is small
    dim.tolphi0 = 1e-6;     % (13/5/24) set phi = 0 when phi <= tolphi
    dim.tolS0 = 0;          % (13/5/24) set S = 0 when S <= tolS0
    dim.tolS1 = 1e-4;       % (13/5/24) set S = 1 when S >= 1-tolS1 (chnaged from 1e-6 to 1e-4 on 09/10/2025)
    
    dim.C = 0.0129*dim.Cstar*dim.P*dim.u;                   % W m^-2 K^-1 (sensible turbulent heat transfer coefficient)
    dim.nu = @(t) (dim.C+4*dim.epsilon*dim.sigma*dim.Tm^3);        % W m^-2 K^-1 (effective turbulent heat transfer coefficient)
    dim.k0 = (dim.dp^2)/180;                                % m^2 (permeability scale)
    dim.s = dim.a/(1-dim.chi*(1-dim.a));                    % dimless (optical property)
    dim.alpha0 = dim.kappa*(1-dim.s)/(1+dim.s);             % m^-1 (CONSTANT absroption coefficient)
    dim.r0 = (dim.kappa^2-dim.alpha0^2)/(2*dim.alpha0);       % m^-1 (CONSTANT scattering coefficient)
    
    dim.Kd = 0; %1;             % dimless (lateral drainage scale, usually the surface slope sin(alpha))
    dim.Ld = 100;           % m (lengthscale over which lateral drainage occurs, eg. lengthscale of catchment until a moulin)
    
    dim.perm = @(phi) dim.k0*phi.^3;                         % m^2 (permeability function k(phi))
    dim.krw = @(S) S.^3;                                    % dimless (relative permeability function krw(S))

    % Microbe and nutrient parameters (22/02/2023) (using mL, NOT L)
    dim.betaA = 5/dim.td; %1/dim.td;          % microbial growth rate (d^-1)
    dim.betaC = 1e-6/dim.td; %1e-5/dim.td;       % rate of consumption of nutrients (micro mol cell^-1 d^-1, 1e-6 in paper)
    dim.alphaPAR = 0.56;    % proportion of solar radiation which is PAR [Jorgensen2001 0.56]
    dim.kPAR = 100;         % saturation radiation (W m^-2) [??]
    dim.kC = 1e-3; %1;             % saturation concentration of nutrients (units m^-3 - > micro mol mL^-1 on 30/07/2024)
    dim.DA = 0;             % microbe diffusion in water
    dim.DC = 0;             % nutrient diffusion in water
    dim.Amax = 10^4;        % maxiumum microbial abundance (cells mL^-1) [Christner2018 10^4]  - updated from 1 units/m^3 on 19/05/2022
    dim.dw = 1e-2/dim.td;          % death rate (water microbes) (d^-1)
    dim.di = 0;             % death rate (ice microbes) (d^-1)

    dim.Ciinf = 1e-3; %1;              % nutrient concentration at depth (micro mol L^-1 -> updated to micro mol mL^-1 on 30/07/2024) - updated from 1 unit/m^3 on 19/05/2022
    dim.Aiinf = 0.01*dim.Amax;  % microbial abundance at depth (cells mL^-1) - updated from 0.5 units/m^3 on 19/05/2022
    
    %% Forcings - most are functions of time (in seconds)
    dim.Qsi = @(t) 200*cos(2*pi*t/dim.td);      % W m^-2 (incoming shortwave radiation)
    dim.Q0 = @(t) -50;                          % W m^-2 (combined longwave radiation and turbulent heat fluxes)
    dim.Rs = @(t) 0;                            % m s^-1 ? (rate of snowfall (solid accumulation))
    dim.Rr = @(t) 0;                            % m s^-1 ? (rate of rainfall)
    dim.Tinf = -10;                        % degrees C (ice temperature at depth - NOT a function of time)
    dim.alpha = @(phi,S,t) dim.alpha0*phi.^0; %*(1+phi.*S);              % m^-1 (absorption coefficient as a function of porosity and saturation)
    dim.r = @(phi,S,t) dim.r0*phi.^0;                                % m^-1 (scattering coefficient as a function of porosity and saturation)
end