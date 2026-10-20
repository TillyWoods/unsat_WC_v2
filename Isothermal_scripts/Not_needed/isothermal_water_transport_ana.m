function [Sana,qana,Zw] = isothermal_water_transport_ana(t,Z,q0,phi,pp)
    % 16th Feb 2024
    % Analytical solution to the isothermal water transport problem WITH
    % lateral drainage, for large Bond number (B>>1, i.e. neglecting
    % surface tension).
    % DIMENSIONLESS.
    %
    % The equation being solved is
    % d/dt(phi*S) + dq/dZ = -qd/Ld = -Kd*q,
    % where q = K*phi^3*S^3 in the B>>1 case.
    %
    % NOTE: if surface tension is included, the RHS is no longer equal to
    % -Kd*q, because the additional surface tension term appears in q but
    % not qd.
    
    K = pp.K;
    Kd = pp.Kd;
    if Kd==0
        Zw = q0^(2/3)*K^(1/3)*t;
    else
        Zw = (3/Kd)*log(1+(1/3)*Kd*K^(1/3)*q0^(2/3)*t);
    end
    Sana = (1./phi).*(q0/K)^(1/3).*exp(-Kd*Z/3).*(Z<Zw);
    qana = q0*exp(-Kd*Z).*(Z<Zw);
end