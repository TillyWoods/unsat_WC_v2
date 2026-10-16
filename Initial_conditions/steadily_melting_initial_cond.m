function ss0 = steadily_melting_initial_cond(Qsi,Q0,pp,dd,options)
    % Calculating the steadily melting solution for Qsi and Q0 
    % (dimensionless), and putting it into an approprate structure to use
    % as an initial condition. By default, the steadily melting solution
    % WITH lateral drainage (as specified by pp.Kd) is used, but there is
    % the option to specify the saturated, no-drainage solution to be used
    % instead, with "drainage"=0.

    arguments
        Qsi (1,1) double    % incoming shortwave radition
        Q0 (1,1) double     % longwave radiation and turbulent heat fluxes
        pp struct           % dimensionless parameters
        dd struct           % discretisation
        options.drainage (1,1) logical = 1   % 1 to for solution with subsurface drainage, 0 for saturated, no-drainage solution
    end

    % The analytical steadily melting state gives us phi,S,T,Zm
    ss0 = analytical_steadily_melting(Qsi,Q0,pp.Tinf,pp,dd);

    if options.drainage == 0
        % Saturated case with no subsurface drainage
        ss0.S = 1*dd.Z.^0;      % saturated
        ss0.q = 0*dd.Zedge.^0;      % no vertical water flux
        ss0.qd = 0*dd.Zedge.^0;     % no subsurface drainage
    else
        % Case with subsurface drainage
        [ss0.S,~,~,~] = sm_drainage_ana_sat(dd.Z,pp,Qsi,Q0);
        [~,ss0.q,ss0.qd,~] = sm_drainage_ana_sat(dd.Zedge,pp,Qsi,Q0);
    end
    
    % Calculate/prescribe the other initial conditions
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp);   % initial enthalpy and total water content
    ss0.h = 0;                          % surface height
    ss0.Zf = 0;                         % position of freezing front (if crust freezes from above)
end