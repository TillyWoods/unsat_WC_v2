function [pw,q,roff,funsat_cell,funsat_edge] = water_fluxes_solve(phi,S,W,hdot,pp,dd,Rs,Rr)
    % Calculates the water fluxes, for a given saturation and porosity
    % solution.

    % 5th June 2024 (modified from water_pressure_and_flux_unsat_flag.m, 13th May 2024)
    % The difference: this new version solves for the water pressure using
    % the same method everywhere, rather than analytically when there is a
    % saturated region above an impermeable base (to allow for the
    % extension to include lateral drainage).
    % FURTHER UPDATE (23rd July 2024): I now solve the equation for q
    % directly in saturated regions bounded on both sides by solid. This
    % works!            

    % Solve for the water pressure pw at the cell centres and the water 
    % fluxes q at the cell edges. A parametrisation is used for pw in the 
    % unsaturated regions and and ODE is solved in the saturated regions.
    % The flux q is calculated differently in satuarated and unsatuarated
    % regions. At the cell edges with one negihbouring cell saturated and
    % one unsaturated, the flux which is largest waway from the saturated
    % cell is used (to help conserve mass).
    %
    % pw = water pressure at cell centres (N column vector)
    %
    % q = water fluxes at cell edges (N+1 column vector)
    %
    % roff = surface runoff (scalar)
    %
    % ss = solution structure containing porosity phi and saturation S at
    % cell centres
    % pp = structure of parameters
    % dd = structure of discretisations (differentiation matrices etc.)
    % Rs = snowfall rate at surface (scalar)
    % Rr = rainfall rate at surface (scalar)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Approximate phi and S on the cell edges, using...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi_edge = harmonic_mean(phi);      % ... harmonic mean for phi
    S_edge = [S(1);S];                  % ... upwinding for phi

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Locate the saturated and unstaurated regions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Tolerances (see usage below)
    tolS1 = pp.tolS1;        % tolerance for saturated region
    tolphi0 = pp.tolphi0;    % tolerance for solid region
    tolS0 = pp.tolS0;        % tolerance for dry region
    
    % Identify regions based on tolerances
    sat = and(S>=1-tolS1,phi>0+tolphi0);    % saturated (and porous)
    dry = and(S<=tolS0,phi>tolphi0);        % dry, i.e. no water (and porous)
    solid = logical(phi<=tolphi0);          % solid (phi approx 0)
    unsat = logical(1-sat-dry-solid);       % unsaturated (and porous)
    unsatd = logical(unsat+dry);     % unsaturated and dry cells
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Identify if a saturated region bounded by solid regions exists
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % If so, the water flux will be solved for directly in these regions,
    % without needing to solve for the water pressure

    sat_sol_edge = edges(sat,solid);    % edges with sat above, solid below
    sol_sat_edge = edges(solid,sat);    % edges with solid above, sat below
    
    % Create an (empty) logical to identify cells that are in a 
    % solid-sat-solid region
    sol_sat_sol = false(size(S));

    % If both a solid-above-sat and sat-above-solid exist then there could 
    % be a solid-sat-solid region
    if and(sum(sat_sol_edge)>0,sum(sol_sat_edge)>0)
        
        % If so, identify...
        sol_sat_ind = find(sol_sat_edge==1);    % indices of solid-above-sat edges
        sat_sol_ind = find(sat_sol_edge==1);    % indices of sat-above-solid edges

        % For each possible bottom of a sol-sat_sol region (i.e. all the
        % sat-above-solid edges), check if there is a continuos set of
        % saturated cells between that and any of the solid-above-sat edges
        % above it. If so, this is a solid-sat-solid region.
        for ii = 1:length(sat_sol_ind)      % for each possible bottom of sol-sat-sol region...  
            ibot = sat_sol_ind(ii);
            itop = find(sol_sat_ind<ibot,1,'first');    % ...find the first possible top of a solo-sat-sol region above it
            if ~isempty(itop)
                if sum((sat(itop:ibot-1)==1))==ibot-itop   % if all the cells between the solid edges are saturated
                    sol_sat_sol(itop:ibot-1)=1;            % set logical identifying solid-sat-solid regions to 1
                end
            end
        end
    end

    % Identify the saturated regions that are NOT bounded by solid on both
    % sides
    sat_not_sol = and(sat,~sol_sat_sol);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve for water pressure pw
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % The water pressure pw is rquired to calculate the water flux q in 
    % (most) saturated regions. As boubdary conditions for calculated pw 
    % in the saturated regions, we also need the water pressure in the 
    % unsaturated regions. 
    
    % pw is calculated differently in unsaturated and
    % saturated regions.

    %%%%%%%%%%%%%%%%%%%%%%%
    % Unsaturated regions
    %%%%%%%%%%%%%%%%%%%%%%%

    % Calculate the water pressure pw in unsaturated regions using 
    % pw = -1/(B*S^2), which we approximate by pw = 0 (since the Bond
    % number B is large and, regardless, including the 1/S^2 term causes
    % numerical issues
    pw = NaN*ones(dd.N,1);
    pw(unsat) = 0;

    % pw is set to 0 in dry regions, i.e. assume they're at atmospheric 
    % pressure
    pw(dry) = 0;  

    % Set pw=0 in solid regions 
    pw(solid) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Saturated regions
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set up linear system to solve 
    % d/dZ(K*phi^3*dpw/dZ)=d/dZ(K*phi^3) + Kd*K*phi^3 
    % for pw in the saturated regions

    k = pp.K*pp.perm(phi_edge);     % permeability (pp.perm(phi) = phi^3 in the equation above)

    % Implement pw=0 at Z=0 when the surface is saturated
    Dcentre = dd.Dcentre;
    Dcentre(1,:) = sparse(1,1,2./dd.dZ(1),1,dd.N);

    % Build linear system
    % (The second term on the rhs is the lateral drainage term)
    Dflux = dd.Dflux;
    lhs = Dflux*spdiags(k,0,dd.N+1,dd.N+1)*Dcentre;
    rhs = Dflux*k + pp.Kd*pp.K*pp.perm(phi);
    
    % Solve in saturated regions NOT BOUNDED BY TWO SOLID REGION 
    % (The term subracted from the rhs ensures continuity of pressure at 
    % interfaces with unsaturated regions)
    pw(sat_not_sol) = lhs(sat_not_sol,sat_not_sol)\(rhs(sat_not_sol)-lhs(sat_not_sol,unsat)*pw(unsat));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate fluxes q at cell edges (between cells of same type)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q = NaN*ones(dd.N+1,1);       % q=NaN when there is no water

    % Flux saturated regions, in terms of pw (and phi)
    qsat = k.*(1-Dcentre*pw); 

    % Flux in unsaturated regions, in terms of S (and phi)
    krw = pp.krw(S_edge);                   % relative permeability
    qunsat = k.*(krw-(2/pp.B)*dd.Dcentre*S);    
    qunsat(1) = Rs+Rr-W(1)*hdot;    % to satisfy (saturated) surface BC

    % If the saturation is below a tolerance, set the flux to zero -
    % NEEDED?
    qunsat(harmonic_mean(S)<pp.tolS0) = 0;
    
    % Identify edges that lie within one region (i.e. same type of cell on
    % both sides)
    solid_edge = edges(solid,solid);                % edge with solid cells on both sides
    sat_not_sol_edge = edges(sat_not_sol,sat_not_sol);             % edge with saturated cells NOT IN A SOL-SAT-SOL REGION on both sides
    unsatd_edge = edges(unsatd,unsatd);             % edge with unsaturated/dry cells on both sides
    unsat_edge = edges(unsat,unsat);                % edges with unsaturated cells on both sides

    % Calculate the flux at these internal cell edges, between 2...
    q(solid_edge) = 0;                      % ...solid cells (no flux)
    q(sat_not_sol_edge) = qsat(sat_not_sol_edge);    % ...saturated cells NOT IN A SOL-SAT-SOL REGION
    q(unsatd_edge) = qunsat(unsatd_edge);   % ...unsaturated/dry cells (use the unsaturated form of q)

    % Solve for q directly (without solving for pw first) in saturated
    % regions bounded by two solid regions, if they exist
    % (This is necessary since the problem for pw in these regions is
    % ill-posed)
    if sum(sol_sat_sol)>0
        bottoms = edges(sol_sat_sol,~sol_sat_sol);  % cells at the bottom of sol-sat-sol regions
        tops = edges(~sol_sat_sol,sol_sat_sol);     % cells at the top of sol-sat-sol regions
        botind = find(bottoms);       % indicies for the bottom cells of sol-sat-sol regions
        topind = find(tops);          % indicies for the top cells of sol-sat-sol regions
        for i = 1:length(botind)    % for each sol-sat-sol region...
            Zbot = dd.Zedge(botind(i)+1);   % edge at bottom
            Ztop = dd.Zedge(topind(i));     % edge at top
            Zspan = [Zbot,Ztop];
            odefunc = @(Z,q) -pp.Kd*interp1(dd.Z,pp.perm(phi),Z);    % dq/dZ = -qd = -Kd*phi^3 in sat region
            q0 = 0;     % no flux at the bottom of the sol-sat-sol region
            [Zout,qout] = ode15s(odefunc,Zspan,q0);
            q(topind(i):botind(i)+1) = interp1(Zout,qout,dd.Zedge(topind(i):botind(i)+1));   % interpolate onto cell edges
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate fluxes q betweeen different types of cell
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Identify the edges between different regions
    sat_unsatd_edge = edges(sat,unsatd);    % edge with saturated above and unsaturated/dry below
    unsatd_sat_edge = edges(unsatd,sat);    % edge with unsaturated/dry above and saturated below
    sol_above = edges(solid,~solid);        % any edge with solid above and not solid below
    sol_below = edges(~solid,solid);        % any edge with solid below and not solid above

    % Calculate the fluxes at the edges between different regions.

    % In order to conserve mass... between saturated and unsaturated/dry 
    % regions, use whichever of the unsaturated (qunsat) and saturated 
    % (qsat) forms of the flux give the largest flux out of the saturated 
    % region
    q(sat_unsatd_edge) = max(qunsat(sat_unsatd_edge),qsat(sat_unsatd_edge));    % flux on edge with saturated above and unsaturated/dry below
    q(unsatd_sat_edge) = min(qunsat(unsatd_sat_edge),qsat(unsatd_sat_edge));    % flux on edge with unsaturated/dry above and saturated below

    q(sol_above) = 0;   % ????
    q(sol_below) = 0;   % no flux at an edge with a solid cell below

    % What do I do about a solid edge above a saturated (or any other) 
    % region? Assume that this will never happen? No, I think it can, but
    % how should I deal with it?

    % Flags to save for the next timestep/iteration:
    funsat_cell = unsat + dry;    % = 1 in unsaturated and dry cells
    funsat_edge = unsat_edge;     % = 1 on edges between to unsaturated cells
    % (EXCLUDES edges neighbouring dry cells - this prevents the porority from going wiggly in dry, porous regions) 

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Calculate surface runoff
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if (3/2)*W(1)-(1/2)*W(2)>=1-pp.tolS1     % if the surface is saturated, we have runoff
       roff = Rs + Rr - q(1) - hdot*W(1);   % (we approximate hdot*W at the surface (cell edge) by the value at the cell centre below)
   else
       roff = 0;
   end

%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%
function edgesab = edges(a,b)
    % Creates a logical identifying the edges where logical a is true on
    % the cell above and logical b is true on the cell below.
    edgesab = and([a(1);a],[b;b(end)]);
end
end
