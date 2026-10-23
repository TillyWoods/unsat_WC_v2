% Make plots from all .mat output files in the folder in which it is run
% 15/10/2026

% Create structure containing details of all files in the directory
listing = dir('**');

for i = 1:length(listing)
    if contains(listing(i).name,'.mat') %,or(contains(listing(i).name,'2026_15'),contains(listing(i).name,'2026_16')))     % only .mat files
        filename = listing(i).name;

        % Load the .mat result file
        load(filename)
        ss=ss_save;

        % % Create names to save plots as
        % base_name = "Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-');

        % Run plotting scripts
        % run('./../../Plotting/melt_runoff_comps.m')
        % run('./../../Plotting/crust_evolution_plot_unsat_periodic_switching.m')
        run('./../../Plotting/internal_surface_periodic_switch_SINGLE.m')
    end
end