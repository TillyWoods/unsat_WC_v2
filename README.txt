Weathering crust code with unsaturated flow and lateral drainage, used in the manuscript "Radiation-driven melting and porous flow in the near-surface layer of an ice sheet" by Tilly Woods and Ian Hewitt.

The details of the model are provided in the manuscript.

To produce the figures in the manuscript, run the following scripts:
	- Fig 2: Isothermal_scripts/FIG2_isothermal_snapshots_3_Kd_vals_fig.m
	- Fig 3: Isothermal_scripts/FIG3_isothermal_ana_num_evol_comp.m
	- Fig 4: Steadily_melting_scripts/FIG4_steadily_Kd_compare_combined_fig.m
	- Fig 6 and 7: Periodic_switching_scripts/FIG6_7_periodic_snapshots_plots.m


To make the plots for figures 5 and 8, data must be manually loaded first:
	1. Load either the drainage or no-drainage solution file:
		- No drainage: Periodic_switching_scripts/wc_periodic_switch_Kd_0_dZ_0-00125_dt_0-0003_14-Oct-2026_18_37_09.mat
		- Drainage: Periodic_switching_scripts/wc_periodic_switch_Kd_20_dZ_0-00125_dt_2-5e-05_23-Oct-2026_22_49_27.mat
	2. Run: ss=ss_save;
	3. Run:
 		- Periodic_switching_scripts/FIG5_crust_evolution_plot_unsat_periodic_switching.m
		- Periodic_switching_scripts/FIG8_melt_runoff_comps

To carry out the simulations described in the manuscript, run the following:
	- RUN_ME_isothermal_run.m
	- RUN_ME_steadily_melting.m
	- RUN_ME_periodic_switching_run.m
Within each of these scripts, it is possible to change the value of the lateral drainage strength Kd.
