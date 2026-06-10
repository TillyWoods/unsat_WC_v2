function ss0 = steadily_melting_initial_cond(Qsi,Q0,pp,dd)
    % Dimensionless

    % The analytical steadily melting state gives us phi,S,T,Zm
    ss0 = analytical_steadily_melting(Qsi,Q0,pp.Tinf,pp,dd);      % use analytical steadily melting state

    [ss0.S,~,~,~] = sm_drainage_ana_sat(dd.Z,pp,Qsi,Q0);
    [~,ss0.q,ss0.qd,~] = sm_drainage_ana_sat(dd.Zedge,pp,Qsi,Q0);
    
    % Calculate/prescribe the other initial conditions
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp);                 % initial enthalpy and total water content
    % ss0.q = 0*dd.Zedge.^0;              % water flux=0
    % ss0.a = pp.a;                       % albedo
    ss0.h = 0;                          % surface height
    ss0.Zf = 0;                         % position of freezing front (if crust freezes from above)

end