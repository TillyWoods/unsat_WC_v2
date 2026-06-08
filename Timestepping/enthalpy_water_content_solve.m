function [Hnew,Wnew,Tnew,phinew,Snew,Tsurfnew,msurfnew,hdotnew] = enthalpy_water_content_solve(ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflag,melt,unsatedge,dd,pp,dt)
    % Build and solve the matrix systems for the enthalpy H and total water
    % content W, with provided forcings and water fluxes. The porosity phi,
    % saturation S and temperature T are recovered, and the melt rate, 
    % temperature and porosity at the surface - msurf, Tsurf, phisurf - 
    % are also calculated, as is the rate of change of surface height hdot.

    % Inputs:
    % ss_old = solution (structure) from previous timestep
    % a = albedo
    % F = subsurface shortwave radiation
    % Rs = snowfall rate
    % Rr = rainfall rate
    % Qsi = incoming shortwave radiation
    % Q0 = longwave radiation and turbulent heat fluxes
    % nu = effective heat transfer coefficient
    % Tinf = surface temperature
    % Hflag = enthalpy flag
    % melt = melt flag
    % unsatedge = unsaturated flag
    % dd = discretisation structure
    % pp = dimensionless parameters structure
    % dt = timestep size

    debug_print = 0;    % Set to 1 to print debugging output
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract values of variables from the previous iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W = ss_old.W;
    phi = ss_old.phi;
    S = ss_old.S;
    q = ss_old.q;
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate permeabilities on cell edges
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Aprroximate phi, S and W on the cell edges using...
    phi_edge = harmonic_mean(phi);  % ...harmonic mean for phi
    S_edge = [S(1);S];              % ...upwinding for S

    k = pp.K*pp.perm(phi_edge);       % permeability from previous iteration (probably phi^3)
    krw = pp.krw(S_edge);        % relative permeability from previous iteration (probably S^3)
    k(1) = 0;           % part of implementing the q = Rs+Rr-W*hdot at Z=0 BC

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Other functions used in the matrix equations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    W_hm = harmonic_mean(W);      % harmonic mean for W
    % (W on its own is approximated on cell edges using the harmonic mean,
    % BUT hdot*W on the cell edges is approximated using upwinding)

    W_inv = 1./W;
    phi_inv = 1./phi;
    phi_inv(or(isinf(phi_inv),isnan(phi_inv))) = 0;
    % (Set 1/phi=0 anywhere it's inf or NaN - all these terms will
    % disappear later anyway since we only need 1/phi in the unsaturated 
    % zones where phi>0)
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create flags for where there is (based on the previous iteration)...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    no_water = (Hflag==1);
    no_ice = (Hflag==3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% H (enthalpy) equation %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start implementing BCs (parts to be done before writing down linear 
    % system) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Surface boundary condition: the boundary condition that we apply
    % depends on the value of the melt flag 'melt'.
    if melt==1      % then the SEB becomes T=0 at Z=0
        dd.Dcentre(1,:) = sparse(1,1,2./dd.dZ(1),1,dd.N);  % using (dT/dZ)_(1/2) approx (T1)/(dZ/2)

    else            % the SEB becomes khat*W*dT/dZ + Upsilon*(chi*(1-a)*Qsi+Q0)-(Upsilon*nuhat+Rs+Rr)*T=0 at Z=0
        dd.Dcentre(1,:) = sparse(1,dd.N);   % since all of (dT/dZ)_(1/2) is now on the RHS
    end

    % Far-field condition (T tends to Tinf)
    dd.Dcentre(dd.N+1,:) = sparse(1,dd.N,-2./dd.dZ(dd.N),1,dd.N);
    dd.cell_below(dd.N+1,:) = sparse(1,dd.N);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the linear system
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M = speye(dd.N) ...
        + dt*dd.Dflux*(min(ss_old.hdot,0)*dd.cell_below ...
        + max(ss_old.hdot,0)*dd.cell_above ...
        + (-pp.khat*spdiags(W_hm,0,dd.N+1,dd.N+1)*dd.Dcentre)*(spdiags((no_water+no_ice).*W_inv,0,dd.N,dd.N))...
        -spdiags(unsatedge.*k,0,dd.N+1,dd.N+1)*(2/pp.B)*dd.Dcentre*spdiags(phi_inv,0,dd.N,dd.N));

    b = ss_old.H ...
        - dt*(dd.Dflux*(pp.khat*spdiags(W_hm,0,dd.N+1,dd.N+1)*dd.Dcentre*(no_ice)*pp.Stefan...
        +pp.Upsilon*F...
        +pp.Stefan*k.*krw.*unsatedge + pp.Stefan*q.*(1-unsatedge))...
        + pp.Stefan*pp.Kd*pp.K.*pp.perm(phi).*pp.krw(S));           % lateral drainage

    % Note that the water flux term q is including implicitly (in terms of
    % H, expressing S = H/(phi*Stefan)) on edges that were between 2
    % unsaturated cells in the previous ITERATION (marked by the unsatedge
    % flag). Elsewhere, the flux q is included explicitly (using the
    % solution from the previous TIMESTEP).
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish implement BCs (parts to be done after writing down linear 
    % system)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Surface boundary condition: the boundary condition that we apply
    % depends on the value of the melt flag 'melt'. The RHS only needs
    % modification in the melt==0 case.

    if melt==0    % the SEB becomes khat*W*dT/dZ + Upsilon*(chi*(1-a)*Qsi+Q0)-(Upsilon*nuhat+Rs+Rr)*T=0 at Z=0  
        b = b - dt*dd.Dflux*sparse(1,1,pp.Upsilon*(pp.chi*(1-a)*Qsi+Q0),dd.N+1,1);      % part of no melt SEB
        M = M - dt*dd.Dflux*(pp.Upsilon*nuhat+Rr+Rs)*sparse([1,1],[1,2],[3/2,-1/2],dd.N+1,dd.N)*(spdiags((no_water+no_ice).*W_inv,0,dd.N,dd.N));            % part of no melt SEB
    end

    % Far-field condition (T tends to Tinf)
    b = b + dt*dd.Dflux*pp.khat*spdiags(W_hm,0,dd.N+1,dd.N+1)*sparse(dd.N+1,1,2*Tinf/dd.dZ(dd.N),dd.N+1,1);
    b = b - dt*dd.Dflux*min(ss_old.hdot,0)*sparse(dd.N+1,1,Tinf,dd.N+1,1);

    % Part of the q = Rs+Rr-W*hdot at Z=0 BC
    % (We approximate hdot*W at the surface (cell edge) by the value at the
    % cell centre below)
    b = b - dt*dd.Dflux*(pp.Stefan*sparse(1,1,Rs+Rr-ss_old.hdot*W(1),dd.N+1,1).*(unsatedge));

    %%%%%%%%%
    % Solve
    %%%%%%%%%
    Hnew = M\b;

    %%%%%%%%%%%%%%%%%
    %% W equation %%
    %%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Set up linear system
    %%%%%%%%%%%%%%%%%%%%%%%%

    M2 = speye(dd.N) + dt*dd.Dflux*(min(ss_old.hdot,0)*dd.cell_below ...
        + max(ss_old.hdot,0)*dd.cell_above);
    
    b2 = ss_old.W - dt*(dd.Dflux*(unsatedge.*k.*(krw-(2/(pp.B*pp.Stefan))*dd.Dcentre*spdiags(phi_inv,0,dd.N,dd.N)*Hnew)...
        + (1-unsatedge).*q)...
        + pp.Kd*pp.K.*pp.perm(phi).*pp.krw(S));             % lateral drainage - DON'T USE k and krw here, use pp.perm, pp.krw


    %%%%%%%%%%%%%%%%%%%%
    % Implement BCs
    %%%%%%%%%%%%%%%%%%%%
    
    % Part of the q = Rs+Rr-W*hdot at Z=0 BC
    % (We approximate hdot*W at the surface (cell edge) by the value at the
    % cell centre below)
    b2 = b2 - dt*dd.Dflux*(sparse(1,1,Rs+Rr-ss_old.hdot*W(1),dd.N+1,1).*(unsatedge));

    % Far-field condition (W tends to 1)
    b2 = b2 - dt*dd.Dflux*min(ss_old.hdot,0)*sparse(dd.N+1,1,1,dd.N+1,1);


    %%%%%%%%%%%
    % Solve
    %%%%%%%%%%%
    Wnew = M2\b2;

    %%%%%%%%%%%%%%%%%%%%%%
    %% Recover T,phi,S
    %%%%%%%%%%%%%%%%%%%%%%
    % Convert from enthalpy and total water content to temperature,
    % porosity and saturation.
    [Tnew,phinew,Snew] = convert_H_to_T(Hnew,Wnew,pp);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate other quantities to output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Surface melt rate, msurf
    msurfnew = (pp.Upsilon*(pp.chi*(1-a)*Qsi+Q0)/pp.Stefan)*melt;
    % (Set to zero if the no melt BC is used)
    
    % Surface temperature, Tsurf
    Tsurfnew = ((3/2)*Tnew(1) - (1/2)*Tnew(2))*(1-melt);
    % (Set to zero if the melt BC is used)
    
    % Surface porosity, phisurf
    phisurf = max((3/2)*phinew(1)-(1/2)*phinew(2),0);
    % (Need the max( ,0) so that it doesn;t become negative when the top 
    % cell is zero and the second is positive)
    
    % Rate of change of surface height, hdot
    hdotnew = (Rs - msurfnew)/(1-phisurf);

    %%%%%%%%%%%%%%
    %% Debugging 
    % %%%%%%%%%%%%

    % Give a warning if the far-field conditions are not satisfied, if the
    % debugging output flag is set to 1
    if debug_print == 1
        if or(or(Wnew(end)<1-1e-5,abs(Hnew(end)-pp.Tinf)>1e-4),unsatedge(end)==1)
            fprintf(['ISSUE!!! Wend = ',num2str(Wnew(end)),', Hend - Tinf = ',num2str(Hnew(end)-pp.Tinf),', unsatedge(end) = ',num2str(unsatedge(end))])
        end
    end

end