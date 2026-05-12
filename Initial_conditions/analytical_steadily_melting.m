function ss = analytical_steadily_melting(Qsi,Q0,Tinf,pp,dd)
    % 27th October 2022 (updated on 8th Nov 2023 to work when there is no 
    % crust but a steadily melting solution exists where the ice is solid 
    % and just the surface is melting (i.e. there is no internal melt)).
    % Dimensionless.

    no_internal_melt = (1-pp.a).*(1-pp.chi).*Qsi./(-Tinf) < ((1-pp.a).*Qsi+Q0)./(pp.Stefan-Tinf);
    
    % Surface melt rate
    ss.msurf = (Tinf.^0)*(pp.Upsilon*(pp.chi*(1-pp.a)*Qsi+Q0)/pp.Stefan.*(1-no_internal_melt) + no_internal_melt*((1-pp.a).*Qsi+Q0)./(pp.Stefan-Tinf));
    
    % Rate of change of surface height
    ss.hdot = -pp.Upsilon*(Qsi*(1-pp.a)+Q0)./(pp.Stefan-Tinf);
    
    % Melting depth
    if no_internal_melt
        ss.Zm = 0;          % Zm=0 if there is no internal melting
    else
        ss.Zm = -(1/pp.kappa)*log(-(Tinf./(pp.Stefan-Tinf)).*((1-pp.a)...
            *Qsi+Q0)./((1-pp.a)*(1-pp.chi)*Qsi));     
    end
    
    % Internal radiation
    ss.F = Qsi*(1-pp.a)*(1-pp.chi).*exp(-pp.kappa*dd.Zedge);
    ss.Fplus = pp.a*Qsi*exp(-pp.kappa*dd.Zedge);
    ss.Fminus = (1-pp.chi*(1-pp.a))*Qsi*exp(-pp.kappa*dd.Zedge);
    
    % Porosity
    ss.phi = (1 + (1/ss.hdot)*(ss.msurf+(pp.Upsilon*Qsi*(1-pp.a)*(1-pp.chi)/pp.Stefan)*(1-exp(-pp.kappa*dd.Z))))...
    .*(dd.Z<ss.Zm); 
    
    % Temperature
    ss.T = (-(pp.Upsilon*Qsi*(1-pp.a)*(1-pp.chi)/(pp.khat*pp.kappa+ss.hdot))*(exp(-pp.kappa*dd.Z)...
        -exp(-pp.kappa*ss.Zm).*exp((ss.hdot/pp.khat)*(dd.Z-ss.Zm))) ...
        + Tinf*(1 - exp((ss.hdot/pp.khat)*(dd.Z-ss.Zm)))).*(dd.Z>ss.Zm);
    
    % Internal melt rate
    ss.mint = -(1/pp.Stefan)*pp.Upsilon*Qsi*(1-pp.a)*(1-pp.chi)*(-pp.kappa*exp(-pp.kappa*dd.Z)).*(dd.Z<ss.Zm);
    
    % Saturation (fully saturated in steadily melting state)
    ss.S = 1*dd.Z.^0;
    
    % Runoff
    ss.roff = -ss.hdot;
end
