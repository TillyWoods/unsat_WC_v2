function dd = discretise_isothermal()
% 24th Feb 2024
% Discretisation to use fo the isothermal water transpoort testing of the
% unsaturated time-dependent weathering crust code.

% Discretisation
Zend = 100;         % size of domain
Zfine = 4;          % size of fine domain near the surface
dZedgeF = 0.01;     % grid size is fine part of domain
dZedgeC = 1e-1;     % grid size in coarse part of domain
NT = 10;            % number of grid cells in which to linearly transition between fine and coarse mesh
Zedge = nonuniform_Zedge(Zend,Zfine,dZedgeF,dZedgeC,NT);

dd = discretise(Zedge);
dd.adaptive = 0;    % no adaptive timestepping
dd.gap = 1000;      % save every 'gap' timesteps
dd.dt = 1e-5;
dd.plot = 0;        % 1 to plot the solution as it's calculated

end