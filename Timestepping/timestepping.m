function ss = timestepping(ss0,pp,dd,tspan,runname)
    % Carry out the simulation, i.e. do the timestepping until the final
    % time is reached (or another condition, e.g. reaching a stead state,
    % is met).

    debug_print = 0;    % Set to 1 to print debugging output

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup for simulation and solution storage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Option to stop simulation if a steady state is reached
    steady_stop = dd.steady_stop;       % 1 to stop once a steady state is reached
    steady_tol = dd.steady_tol;         % tolerance used to determine steady state

    % Option to stop simulation if periodic solution is reached
    periodic_stop = dd.periodic_stop;   % the size of the period to stop once a periodic state is reached, 0 if not
    periodic_tol = dd.periodic_tol;     % tolerance used to determine periodic state (percentage error (old-new)/new): 0.01 means 1% error)

    % Decide how often to output the solution (i.e. save as .mat)
    save_period = 60;           % output every save_period time
    savetime = save_period;     % will be updated every output to give the time for the next output
    part = 0;                   % counter for number part of solution we're saving as we go

    % Extract time span and itialise time
    tstart = tspan(1);      % start time
    tend = tspan(2);        % end time
    t = tstart;

    j = 0;              % timestep counter
    
    %%%%%%%%%%%%%%%%%%%%%%
    %% Initial condition
    %%%%%%%%%%%%%%%%%%%%%%
    % Some of the variables that will be outputted later are not provided
    % by the inputted initial condition so are presribed here.
    ss_step = ss0;          % inputted initial condition
    ss_step.crust = ss_step.Zm-ss_step.Zf>1e-12;    % logical identifying crust location
    ss_step.t = tstart;                     % time
    ss_step.Z = dd.Z;                       % depth below surface (cell centres)
    ss_step.Zedge = dd.Zedge;               % depth below surface (cell edges)
    ss_step.Tsurf = (3/2)*ss_step.T(1) - (1/2)*ss_step.T(2);    % surface temperature - linear extrapolation
    ss_step.pw = NaN*dd.Z.^0;               % water pressure
    ss_step.qd = pp.K*pp.Kd*pp.perm(ss_step.phi).*pp.krw(ss_step.S);   % lateral drainage
    ss_step.j = j;                          % timestep counter
    ss_step.i = 1;                          % saved solution counter

    % Prescribed forcings to be outputted
    ss_step.Qsi = pp.Qsi(t);        % incoming shortwave radiation
    ss_step.Q0 = pp.Q0(t);          % longwave radiation and turbulent heat fluxes
    ss_step.Tinf = pp.Tinf;         % temperature at depth
    ss_step.Rs = pp.Rs(t);          % snowfall rate
    ss_step.Rr = pp.Rr(t);          % rainfall rate
    ss_step.nu = pp.nuhat(t);       % effective heat transfer coefficient

    % Initialise cumulative quantities to zero
    ss_step.msurf_cum = 0;              % cumulative surface melt
    ss_step.mint_cum = 0*dd.Z.^0;       % cumulative internal melt
    ss_step.freeze_cum = 0*dd.Z.^0;     % cumulative internal freezing
    ss_step.roff_cum = 0;               % cumulative surface runoff
    ss_step.qd_cum = zeros(size(dd.Z));    % cumulative lateral drainage

    % Initialise flags (unsat flag, melt flag, enthalpy flag)
    S = ss_step.S;
    phi = ss_step.phi;
    ss_step.unsatedge = and(and([S(1);S]<=1-pp.tolS1,...
        [phi(1);phi]>=pp.tolphi0),and([S;S(end)]<=1-pp.tolS1,...
        [phi;phi(end)]>=pp.tolphi0));       % unsat flag
    ss_step.melt = ss_step.msurf>0;                     % melt flag
    ss_step.Hflag = H_flag(ss_step.H,ss_step.W,pp);     % enthalpy flag

    % Store the initial condition
    ss(1) = ss_step;

    %%%%%%%%%%%%%%%%
    %% Timestepping
    %%%%%%%%%%%%%%%%
    warning('No warning.')
    while and(t<tend,1)
        % Set the timestep size using the constant value stored in dd
        if j<1000
            dt = 0.025*dd.dt;
        else
            dt = dd.dt;
        end
        dt = min(dt,tend-t);    % adjust final timestep to hit end time
        t = t + dt;             % time at new timestep
        j = j + 1;              % timestep number
        
        % Evaluate the forcings at the current timestep
        Qsi = pp.Qsi(t);        % incoming shortwave radiation
        Q0 = pp.Q0(t);          % longwave radiation adn turbulent heat fluxes
        Rs = pp.Rs(t);          % snowfall rate
        Rr = pp.Rr(t);          % rainfall rate
        nu = pp.nuhat(t);       % effective heat transfer coefficient
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Carry out timestepping
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        ss_step = SEB_iterations(ss_step,Qsi,Q0,pp.Tinf,Rs,Rr,nu,pp,dd,dt);
        if debug_print == 1
            fprintf([newline,'t = ',num2str(t),newline])
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store the timesteps we want to store
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if or(mod(j-1,dd.gap)==0,t>=tend)   % tend bit makes sure we always store the final timestep
            if debug_print == 0
                fprintf([newline,'t = ',num2str(t),newline])
            end
            
            % Update index i of stored solution steps
            if t>=tend
                i=i+1;
            else
            i = (j-1)/dd.gap + 2;
            end
            ss_step.i = i;                % stored solution number
            ss_step.t = t;                % time
            ss_step.j = j;                % timestep number
            ss_step.Z = dd.Z;             % positions of cell centres
            ss_step.Zedge = dd.Zedge;     % positions of cell edges
            ss(i) = ss_step;

            % Set phi and S to 0 or 1 below or above the tolerances in the
            % solution to save (but not in the solution to use in the next
            % timestep)
            ss(i).phi(ss(i).phi<=tolphi0) = 0;      % solid
            ss(i).S(ss(i).S<=tolS0) = 0;            % dry
            ss(i).S(ss(i).S>=1-tolS1) = 1;          % saturated

            % Set S to 0 in solid regions for clarity in plots
            ss(i).S(ss(i).phi<=pp.tolphi0) = 0;

            % Save solution structure every 'save_period' time (if zero, 
            % only save at the end of the run)
            if save_period > 0
                if t > savetime
                    part = part + 1;            % number part of solution that we're saving
                    save_time = 120; %inf;      % save the last 'save_time' amount of time of data
                    save_ind = ss(end).t-[ss.t]<save_time;
                    ss_save = ss(save_ind);

                    save(["Outputs/wc_periodic_switch_"+runname+"_"+strrep(strrep(string(date),' ','_'),':','_')+".mat"],"ss_save","pp","dd")
                    
                    % Add another period to savetime so that we don't save
                    % any more solutions until another period has passed
                    savetime = savetime + save_period;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot solution (if desired, e.g. for debugging)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if dd.plot==1
                figure(1)
                clf;
                width = 2;
                hold on
                p1 = plot(ss(i).phi,dd.Z,'k-','linewidth',width);
                p2 = plot(ss(i).S,dd.Z,'b-','linewidth',width);
                p3 = plot(10*ss(i).q,dd.Zedge,'g-','linewidth',width);
                p4 = plot(ss(i).pw,dd.Z,'r-','linewidth',width);
                p5 = plot(ss(i).T,dd.Z,'m-','linewidth',width);
                p6 = plot(ss(i).F,dd.Zedge,'c-','linewidth',width);
                p9 = plot(ss(i).qd,dd.Z,'color',[0.4,0.2,0.5],'linewidth',width);
                p7 = plot([0,1],[ss(i).Zm,ss(i).Zm],'k--','linewidth',width);
                % p8 = plot(ss(i).H,dd.Z,'color',0.5*[1 1 1],'linewidth',2*width,'marker','x');
                set(gca,'YDir','reverse')
                legend([p1,p2,p3,p4,p5,p6,p9,p7],{'$\phi$','$S$','$10q$','$p_{w}$','$T$','$F$','$q_\mathrm{d}$','$Z_\mathrm{m}$'},'interpreter','latex','location','southeast');
                title(['$t=$',num2str(t,'%.2f')],'interpreter','latex');
                ylabel('$Z$','interpreter','latex');
                xlim([-0.2,1.2]) %xlim([-0.5,0.1]);
                ylim([0,0.5]) %ylim([0,3])
                hold off
                shg
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Conditions to stop the timestepping early:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % If MATLAB gives a singular or badly scaled matrix warning, stop
        % the timestepping:
        [~,msgid] = lastwarn;
        if or(strcmp('MATLAB:nearlySingularMatrix',msgid),strcmp('MATLAB:SingularMatrix',msgid))
            disp('Timestepping stopped due to singular or badly scaled matrix.')
            break
        end
        
        % If desired, stop timestepping when steady state is reached
        if and(steady_stop,norm(ss(i).T-ss(i-1).T)+norm(ss(i).phi-ss(i-1).phi)+norm(ss(i).S-ss(i-1).S)<steady_tol)
            disp('Steady state reached.')
            break
        end
        
        % If desired, stop timestepping when perioic state is reached
        np = ceil(periodic_stop/dd.dt/dd.gap);
        two_periods = and(length([ss.t])>2*np,mod(length([ss.t]),np)==1);      % if an integer number of at least 2 periods have happened
        if two_periods
            norm_diff = norm([ss(i-np:i).T]-[ss(i-2*np:i-np).T])+norm([ss(i-np:i).phi]-[ss(i-2*np:i-np).phi])...
                +norm([ss(i-np:i).S]-[ss(i-2*np:i-np).S]);      % norm of the difference between the latest two periods of solution
            norm_one = norm([ss(i-np:i).T])+norm([ss(i-np:i).phi])...
                +norm([ss(i-np:i).S]);      % norm of the latest period of solution
            errorp = norm_diff/norm_one    % percentage error
            small_diff = errorp < periodic_tol;      % if the difference between the last 2 periods is small
        else
            small_diff=0;
        end
        if and(periodic_stop>0,and(two_periods,small_diff))
            break
        end
    end 
end
