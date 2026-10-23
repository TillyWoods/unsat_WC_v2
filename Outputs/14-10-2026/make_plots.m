% Make plots from all .mat output files in the folder in which it is run
% 15/10/2026

% Create structure containing details of all files in the directory
listing = dir('**');

for i = 1:length(listing)
    if and(contains(listing(i).name,'.mat'),contains(listing(i).name,'Kd_20'))     % only .mat files
        filename = listing(i).name;

        % Load the .mat result file
        load(filename)
        ss=ss_save;

        % % Create names to save plots as
        % base_name = "Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-');

        % Run plotting scripts
        % run('./../../Plotting/melt_runoff_comps.m')
        % run('./../../Plotting/crust_evolution_plot_unsat_periodic_switching.m')
        % run('./../../Plotting/internal_surface_periodic_switch_SINGLE.m')

        plot_times_freeze = [22,22.4,22.8,23,24,25];
        plot_times_melt = [42,43,43.2,43.5,43.7,44]; %[44,45,45.5,46,48,51];
        output_dir = "./../../Plotting/OutputPlots_14-10-2026";
        %plotting_snapshots_periodic_switching_func(ss,pp,plot_times_freeze,output_dir=output_dir)
        plotting_snapshots_periodic_switching_func(ss,pp,plot_times_melt,output_dir=output_dir)
    end
end