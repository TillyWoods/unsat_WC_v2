function [H,W,T,phi,S,Tsurf,msurf,hdot] = one_iteration_HW_new_iter_no_F_it(ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflag,melt,unsatedge,dd,pp,dt)
    % [H,W,F,Fplus,Fminus,a,Tsurf,msurf,unsatedge]

    % 26-08-2024 - modified to calculate F once at the start of each timestep,
    % rather than inside iterations. Only the enthalpy and unsat flags are
    % iterated.

    % F and a are calculated at the start of the timestep.

    % ss_old = solution from previous timestep
    
    W = ss_old.W; %W_it; %ssit.W;
    phi = ss_old.phi; %phi_it;
    S = ss_old.S; %S_it;
    % F = F_it; %ssit.F;
    q = ss_old.q; % q_it; %ssit.q;
    % a = ra_it;
    % qd = qd_it;
    
    % Harmonic means of phi, S, W from the orevious iteration
    phi_hm = harmonic_mean(phi);
    S_hm = harmonic_mean(S);
    W_hm = harmonic_mean(W);

    phi_edge = phi_hm;
    S_edge = [S(1);S];      % upwinding

    W_inv = 1./W;

    phi_inv = 1./phi;
    phi_inv(or(isinf(phi_inv),isnan(phi_inv))) = 0;   % set 1/phi=0 anywhere it's inf or NaN - all these terms will disappear later anyway since we only need 1/phi in the unsaturated zones where phi>0
    
    % Flux in unsaturated region
    k = pp.K*phi_edge.^3;                             % permeability from previous iteration
    krw = S_edge.^3; %S_hm.^3;                                  % relative permeability from previous iteration
    k(1) = 0;           % part of the q = Rs+Rr-W*hdot at Z=0 BC (added 15/05/2024)
    
    no_water = (Hflag==1);
    no_ice = (Hflag==3);

    %% H equation

    % Implement BCs (parts to be done before writing down linear system) 
    if melt==1      % SEB becomes T=0 at Z=0
        dd.Dcentre(1,:) = sparse(1,1,2./dd.dZ(1),1,dd.N);  % since Tsurf = 0
        % If the code was extended to allow for flow in pure water
        % regions, the first row of the q.cell_above matrix would need
        % to be set to zero too.
    else            % SEB becomes khat*W*dT/dZ + Upsilon*(chi*(1-a)*Qsi+Q0)-(Upsilon*nuhat+Rs+Rr)*T=0 at Z=0
        dd.Dcentre(1,:) = sparse(1,dd.N);
    end
    dd.Dcentre(dd.N+1,:) = sparse(1,dd.N,-2./dd.dZ(dd.N),1,dd.N);   % part of the T tends to Tinf BC   

    dd.cell_below(dd.N+1,:) = sparse(1,dd.N);    % part of the T tends to Tinf BC (ADDED 09/06/2023)

    % Whm = harmonic mean of W from previous iteration

    % Set up the linear system
    M = speye(dd.N) + dt*dd.Dflux*(min(ss_old.hdot,0)*dd.cell_below ...
        + max(ss_old.hdot,0)*dd.cell_above + (-pp.khat*spdiags(W_hm,0,dd.N+1,dd.N+1)*dd.Dcentre)*(spdiags((no_water+no_ice).*W_inv,0,dd.N,dd.N))...
        -spdiags(unsatedge.*k,0,dd.N+1,dd.N+1)*(2/pp.B)*dd.Dcentre*spdiags(phi_inv,0,dd.N,dd.N));
    b = ss_old.H - dt*(dd.Dflux*(pp.khat*spdiags(W_hm,0,dd.N+1,dd.N+1)*dd.Dcentre*(no_ice)*pp.Stefan...
        +pp.Upsilon*F...
        +pp.Stefan*k.*krw.*unsatedge + pp.Stefan*q.*(1-unsatedge))...
        + pp.Stefan*pp.Kd*pp.K.*pp.perm(phi).*pp.krw(S)); % pp.Stefan*qd); %qd(1:end-1)); %             % lateral drainage - DON'T USE k and krw here, use pp.perm, pp.krw
        % + pp.Stefan*pp.Kd*pp.K.*pp.perm(phi_edge(1:end-1)).*pp.krw(S_edge(1:end-1)));

    % Implement BCs (parts to be done after writing down linear system) 
    if melt==0
%             M = M + dt*dd.Dflux*sparse(1,1,pp.Upsilon*(pp.chi*(1-a)*Qsi+Q0),dd.N+1,1);
        b = b - dt*dd.Dflux*sparse(1,1,pp.Upsilon*(pp.chi*(1-a)*Qsi+Q0),dd.N+1,1);      % part of no melt SEB
        M = M - dt*dd.Dflux*(pp.Upsilon*nuhat+Rr+Rs)*sparse([1,1],[1,2],[3/2,-1/2],dd.N+1,dd.N)*(spdiags((no_water+no_ice).*W_inv,0,dd.N,dd.N));            % part of no melt SEB
    end
    b = b + dt*dd.Dflux*pp.khat*spdiags(W_hm,0,dd.N+1,dd.N+1)*sparse(dd.N+1,1,2*Tinf/dd.dZ(dd.N),dd.N+1,1);   % part of the T tends to Tinf BC
    b = b - dt*dd.Dflux*min(ss_old.hdot,0)*sparse(dd.N+1,1,Tinf,dd.N+1,1);      % part of T tends to Tinf BC (ADDED 09/06/2023)
    b = b - dt*dd.Dflux*(pp.Stefan*sparse(1,1,Rs+Rr-ss_old.hdot*W(1),dd.N+1,1).*(unsatedge));      % part of the q = Rs+Rr-W*hdot at Z=0 BC (we approximate hdot*W at the surface (cell edge) by the value at the cell centre below) (added 15/05/2024)
    % b = b + dt*dd.Dflux*(-pp.khat*spdiags(W_hm,0,dd.N+1,dd.N+1)*dd.Dcentre)*((no_ice)*pp.Stefan);    % 30/11/2023 - to allow pure water to exceed the melting temp

    % Solve
    Hnew = M\b;

    %% W equation
    % Set up linear system
    M2 = speye(dd.N) + dt*dd.Dflux*(min(ss_old.hdot,0)*dd.cell_below ...
        + max(ss_old.hdot,0)*dd.cell_above);
    b2 = ss_old.W - dt*(dd.Dflux*(unsatedge.*k.*(krw-(2/(pp.B*pp.Stefan))*dd.Dcentre*spdiags(phi_inv,0,dd.N,dd.N)*Hnew)...
        + (1-unsatedge).*q)...
        + pp.Kd*pp.K.*pp.perm(phi).*pp.krw(S)); %qd); %qd(1:end-1)); %             % lateral drainage - DON'T USE k and krw here, use pp.perm, pp.krw
        % + pp.Kd*pp.K.*pp.perm(phi_edge(1:end-1)).*pp.krw(S_edge(1:end-1)));

    % Implement BCs
    b2 = b2 - dt*dd.Dflux*(sparse(1,1,Rs+Rr-ss_old.hdot*W(1),dd.N+1,1).*(unsatedge));    % part of the q = Rs+Rr-W*hdot at Z=0 BC (we approximate hdot*W at the surface (cell edge) by the value at the cell centre below) (added 15/05/2024)
    b2 = b2 - dt*dd.Dflux*min(ss_old.hdot,0)*sparse(dd.N+1,1,1,dd.N+1,1);      % part of W tends to 1 BC (ADDED 20/05/2024)

    % Solve
    Wnew = M2\b2;

    % Wnew(Wnew>=1-pp.tolphi0) = 1;     % don't do this! Causes things not to be conserved

    %% Recover T,phi,S
    [Tnew,phinew,Snew] = convert_H_to_T(Hnew,Wnew,pp);
    if sum(Snew(dd.Z<2)>=1)>=1
        Snew;
    end

    if or(or(Wnew(end)<1-1e-5,abs(Hnew(end)-pp.Tinf)>1e-4),unsatedge(end)==1)
        fprintf(['ISSUE!!! Wend = ',num2str(Wnew(end)),', Hend - Tinf = ',num2str(Hnew(end)-pp.Tinf),', unsatedge(end) = ',num2str(unsatedge(end))])
    end

    %% Calculate Tsurf, msurf, phisurf, hdot
    msurfnew = (pp.Upsilon*(pp.chi*(1-a)*Qsi+Q0)/pp.Stefan)*melt;  % msurf is set to zero if the no melt BC is used.
    Tsurfnew = ((3/2)*Tnew(1) - (1/2)*Tnew(2))*(1-melt);                 % Tsurf is set to zero if the melt BC is used
    phisurf = max((3/2)*phinew(1)-(1/2)*phinew(2),0);               % surface porosity - need the max( ,0) so that it doesn;t become negative when the top cell is zero and the second is positive
    hdotnew = (Rs - msurfnew)/(1-phisurf);  % rate of change of surface height

    H = Hnew;
    W = Wnew;
    % a = anew;
    Tsurf = Tsurfnew;
    msurf = msurfnew;
    hdot = hdotnew;
    T = Tnew;
    phi = phinew;
    S = Snew;
end
