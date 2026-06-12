function dd = discretise_nonuniform()
% 3-layered nonuniform grid discretisation with a very fine layer near the
% surface, a fine layer below, and coarse below that, with the grid spacing
% being uniform within each of those regions.

% Discretisation
Zend = 100;         % size of domain
Zvfine = 1.5;       % size of very fine domain near the surface
dZedgevF = 0.0025;  % grid size in very fine part of domain
Zfine = 4;          % bottom of the fine domain that extends below the very fine
dZedgeF = 0.01;     % grid size is fine part of domain
dZedgeC = 1e-1;     % grid size in coarse part of domain
NT = 10;            % number of grid cells in which to linearly transition between fine and coarse mesh
Zedge = nonuniform_Zedge_3_parts(Zend,Zvfine,dZedgevF,Zfine,dZedgeF,dZedgeC,NT);

dd = discretise(Zedge);
dd.adaptive = 0;    % no adaptive timestepping
dd.gap = 800;      % save every 'gap' timesteps
dd.dt = 0.5*0.25*1e-3;
dd.plot = 0;        % 1 to plot the solution as it's calculated

end