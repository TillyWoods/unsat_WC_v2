Semi-implicit water flux weathering crust code with unsaturated flow and lateral drainage.

The code in this folder solves the above problem with:
- the water fluxes between two unsaturated cells included implicitly

- iterations carried out to update the enthalpy and unsaturated flags (with all other non-linear variables included explicitly, i.e. from the previous timestep, not the previous iteration)

- internal radiation F is calculated analytically at each timestep using a scalar extinction coefficient

- F is calculated outside of the iterations, at the start of each timestep (in theory using phi and S from the previous timestep if the code was generalised to allow radiative properties to vary with depth)

- the water pressure is solved for in saturated regions, using pw=0 as the boundary condition with the surface AND with unsaturated regions

- in saturated regions with two solid boundaries, q is solved for directly using ode15s - I should maybe change this to solve it manually and to apply to any region with a solid base

- the melt flag is updated separately from the enthalpy and unsaturated flag iterations - the melt flag from the previous timestep is assumed, one set of iterations is carried out, then if the melt flag could have be wrong then the set of iterations is redone using the other melt flag
