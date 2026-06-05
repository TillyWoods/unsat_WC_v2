function ss = timestepping_unsat_implicit(ss0,pp,dd,tspan)
    % 27th October 2022

    debug_print = 0;    % Set to 1 to print debugging output
    
    stopZm = dd.stopZm;
    ind_stop = stopZm; %or(isstring(tspan),ischar(tspan));

    steady_stop = dd.steady_stop;       % 1 to stop once a steady state is reached
    steady_tol = dd.steady_tol;         % tolerance used to determine steady state

    periodic_stop = dd.periodic_stop;   % the size of the period to stop once a periodic state is reached, 0 if not
    periodic_tol = dd.periodic_tol;     % tolerance used to determine periodic state (percentage error (old-new)/new): 0.01 means 1% error)

    % if ind == 1
    %     tstart = 0;
    %     tend = 1e10;
    % else
    %     tstart = tspan(1);      % start time
    %     tend = tspan(2);        % end time
    % end

    tstart = tspan(1);      % start time
    tend = tspan(2);        % end time

    % if isstring(tend)==0
    %     stop_cond = (ss_step.Zm<1e-6);
    % else
    %     stop_cond = (t<tend);
    % end
  
    t = tstart;
    ss_step = ss0;      % initial condition
    ss_step.crust = ss_step.Zm-ss_step.Zf>1e-12;
    ss_step.melt = ss_step.msurf>0; %%%%%%
    ss_step.t = tstart;
    ss_step.Z = dd.Z;
    ss_step.Zedge = dd.Zedge;
    ss_step.Tsurf = (3/2)*ss_step.T(1) - (1/2)*ss_step.T(2);
    ss_step.pw = NaN*dd.Z.^0;
    ss_step.j = 0;
    ss_step.i = 1;

    ss_step.Qsi = pp.Qsi(t);
    ss_step.Q0 = pp.Q0(t);
    ss_step.Tinf = pp.Tinf;
    ss_step.Rs = pp.Rs(t);
    ss_step.Rr = pp.Rr(t);
    ss_step.alpha = pp.alpha(ss_step.phi,ss_step.S,t);
    ss_step.r = pp.r(ss_step.phi,ss_step.S,t);
    ss_step.nu = pp.nuhat(t);

    ss_step.msurf_cum = 0;      % cumulative surface melt (ADDED 29/01/2024)
    ss_step.mint_cum = 0*dd.Z.^0;       % cumulative internal melt (ADDED 29/01/2024)
    ss_step.freeze_cum = 0*dd.Z.^0;     % cumulative internal freezing (ADDED 09/08/2024)
    ss_step.roff_cum = 0;               % cumulative surface runoff (ADDED 17/09/2024)

    ss_step.qd = pp.K*pp.Kd*pp.perm(ss_step.phi).*pp.krw(ss_step.S);   % lateral drainage (ADDED 13/03/2024, with lengthscale over which drainage occurs absorbed into the dimensionless form of qd)
    ss_step.qd_cum = zeros(size(dd.Z));         % cumulative lateral drainage (ADDED 13/03/2024)

    ss_step.Hflag = H_flag(ss_step.H,ss_step.W,pp);
    S = ss_step.S;
    phi = ss_step.phi;
    ss_step.unsatedge = and(and([S(1);S]<=1-pp.tolS1,[phi(1);phi]>=pp.tolphi0),and([S;S(end)]<=1-pp.tolS1,[phi;phi(end)]>=pp.tolphi0));
    ss_step.melt = ss_step.msurf>0;

    ss(1) = ss_step;        % store the initial condition
    
    j = 0;
    crust_gone = 0;
    
    warning('No warning.')
    while and(t<tend,1) %logical(1-ind_stop*(crust_gone>=dd.gap)))              % if stopZm=1, timestepping stops when Zm<=1-12 for dd.gap timesteps in a row
    % while and(t<tend,logical(1-ind*((ss_step.Zm-ss_step.Zf)<=1e-12))) % if stopZm=1, timestepping stops when Zm<=1e-12
    % while (ind*(ss_step.Zm>1e-6)+(1-ind)*(t<tend))==1 %t<tend       % if stopZm=1, timestepping stops when Zm<=1e-6
        % Set the timestep size using the constant timestep size stored in
        % dd.
        dt = dd.dt; %dt_adaptive(ss_step,pp,dd)*(dd.adaptive==1) + dd.dt*(1-dd.adaptive==1);
        dt = min(dt,tend-t);
        t = t + dt; %  t = min(t + dt,tend);  %      % time at new timestep
        j = j + 1;          % timestep number
        
        % Evaluate the forcings at the current timestep
        Qsi = pp.Qsi(t);
        Q0 = pp.Q0(t);
        Rs = pp.Rs(t);
        Rr = pp.Rr(t);
        alpha = @(phi,S) pp.alpha(phi,S,t);
        r = @(phi,S) pp.r(phi,S,t);
        nu = pp.nuhat(t);
        
        % Carry out timestepping
        ss_step = full_timestep_implicit_qflux(ss_step,Qsi,Q0,pp.Tinf,Rs,Rr,alpha,r,nu,pp,dd,dt);
        if debug_print == 1
            fprintf([newline,'t = ',num2str(t),newline])
        end

        crust_gone = double(crust_gone) + double((ss_step.Zm-ss_step.Zf)<=1e-12);

        % Store the timesteps we want to store
        if or(mod(j-1,dd.gap)==0,t>=tend)   % tend bit makes sure we always store the final timestep
            if debug_print == 0     % print time only at timesteps to save
                fprintf([newline,'t = ',num2str(t),newline])
            end
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

            % Plot solution
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
        

        % if (crust_gone>=dd.gap)==1
        %     ind_stop
        % end
    end 
end