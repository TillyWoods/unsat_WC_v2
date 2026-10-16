function dd = discretise_isothermal(options)
% Discretisation to use for the isothermal water transport test case (and 
% the steadily melting test case) of the unsaturated weathering crust code
% with lateral drainage.

arguments
    options.dt (1,1) double = 0.25*1e-2;    % timestep size
end

% Discretisation
Zend = 100;         % size of domain
Zfine = 4;          % size of fine domain near the surface
dZedgeF = 0.01;     % grid size is fine part of domain
dZedgeC = 1e-1;     % grid size in coarse part of domain
NT = 10;            % number of grid cells in which to linearly transition between fine and coarse mesh
Zedge = nonuniform_Zedge(Zend,Zfine,dZedgeF,dZedgeC,NT);

dd = discretise(Zedge);
dd.adaptive = 0;    % no adaptive timestepping
dd.dt = options.dt; % timestep size
dd.gap = floor(0.25*0.1/dd.dt);      % save every 'gap' timesteps
dd.plot = 0;        % 1 to plot the solution as it's calculated

end