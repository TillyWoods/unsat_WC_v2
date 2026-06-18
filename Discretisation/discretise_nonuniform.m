function dd = discretise_nonuniform(options)
% 3-layered nonuniform grid discretisation with a very fine layer near the
% surface, a fine layer below, and coarse below that, with the grid spacing
% being uniform within each of those regions.

arguments
% Discretisation (with default values)
options.Zend (1,1) double = 100;         % size of domain
options.Zvfine (1,1) double = 1.5;       % size of very fine domain near the surface
options.dZedgevF (1,1) double = 0.005;   % grid size in very fine part of domain
options.Zfine (1,1) double = 4;          % bottom of the fine domain that extends below the very fine
options.dZedgeF (1,1) double = 0.01;     % grid size is fine part of domain
options.dZedgeC (1,1) double = 1e-1;     % grid size in coarse part of domain
options.NT (1,1) double = 10;            % number of grid cells in which to linearly transition between fine and coarse mesh

options.dt (1,1) double = 0.5*0.25*1e-3;    % timestep size
end

Zedge = nonuniform_Zedge_3_parts(options.Zend,options.Zvfine,...
    options.dZedgevF,options.Zfine,options.dZedgeF,options.dZedgeC,options.NT);

dd = discretise(Zedge);
dd.adaptive = 0;    % no adaptive timestepping
dd.dt = options.dt;
dd.gap = 400/(dd.dt/(0.25*1e-3));      % save every 'gap' timesteps
dd.plot = 0;        % 1 to plot the solution as it's calculated

end