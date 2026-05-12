function dd = discretise(Zedge)
    % 19th Oct 2022
    % Zedge is a column vector of the of length N+1 containing the
    % positions of the cell ends (from 0 to 'infinity'). From this, all the
    % other discretisations we might want are made. The timestep size dt is
    % also astored here.
    
    % Timestepping
    dd.dt = 1e-3;                       % timestep
    dd.adaptive = 1;                    % 1 to use adaptive timestepping, 0 to use dd.dt every timestep
    dd.gap = 10;                        % store every 'gap' timesteps
    dd.save_mean = 0;                   % 1 to save the mean solution over 'gap' timesteps, rather than the snapshot at the 'gap'th timestep
    dd.stopZm = 0;                      % = 1 to stop the simulation when the crust thicknes becomes zero

    dd.steady_stop = 0;                 % 1 to stop when a steady state is reached (norm(T-Told)+norm(S-Sold)+norm(phi-phiold))<steady_tol
    dd.steady_tol = 1e-4;

    dd.periodic_stop = 0;               % the size of the period to stop once a periodic state is reached, 0 if not
    dd.periodic_tol = 1e-3;
    
    % Grid
    dd.Zedge = Zedge(:);                % positions of cell edges (length N+1) (the colon ensures that dd.Zedge is a column vector)
    dd.N = length(dd.Zedge)-1;          % number of cells (one less than number of edges)
    dd.dZ = diff(dd.Zedge);             % cell widths (length N)
    dd.Z = dd.Zedge(1:end-1)+dd.dZ/2;   % positions of cell centres (length N)
    dd.dZcentre = diff(dd.Z);           % distances between cell centres (length N-1)
    
    % Differentiation matrices
    dd.Dflux = spdiags([-1./dd.dZ,1./dd.dZ],[0,1],dd.N,dd.N+1); % acts on cell edges to give result on cell centres
    dd.Dcentre = spdiags([-1./dd.dZcentre;0],-1,dd.N+1,dd.N) + spdiags([0;1./dd.dZcentre],0,dd.N+1,dd.N);
    
    % Other matrices
    dd.cell_below = spdiags(ones(dd.N,1),0,dd.N+1,dd.N);     % approximates a cell edge value by the value of the cell below (size (N+1)xN)
    dd.cell_above = spdiags(ones(dd.N,1),-1,dd.N+1,dd.N);    % approximates a cell edge value by the value of the cell above (size (N+1)xN) 
end