function [H,W,T,phi,S,F,a,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag,melt] = SEB_iterations(ss_old,Rs,Rr,Qsi,Q0,nuhat,Tinf,dd,pp,dt)
% Iterate the full enthalpy, water content annd water fluxes solutions
% until the correct (i.e. consistent) surface energy balance is reached.
% (If the assumed SEB works (unless msurf = 0 or Tsurf = 0), stick with it.
% If not (i.e. if we get negative or zero surface melt, or positive or zero
% surface temperature), try the other case. If neither work (hopefully 
% never the case) or both work, keep the second attempt.)


debug_print = 0;    % Set to 1 to print debugging output

meltold = ss_old.melt;
melt = meltold;
Hflagold = ss_old.Hflag;
unsatedgeold = ss_old.unsatedge;

%% Calculate internal shortwave radiation F
% This is done once each timestep, for the given value of incoming
% shortwave radiation Qsi at that timestep.
% (This could be modified to allow for radiative transfer properties that
% also evolve over time).

kappa = pp.kappa;       % extinction coefficient
a = pp.a;               % albedo
chi = pp.chi;           % shortwave partitioning parameter

F = (1-chi)*(1-a)*Qsi*exp(-kappa*dd.Zedge);

 %% Step 1: carry out the timestep...
    % ... assuming that the melt and enthlapy flags are that of the
    % previous timestep.
    % NOTE: within this, the a number (probably 10 - see function below to 
    % be sure) of iterations are carried out in which the enthalpy flag is 
    % updated each time if it was incorrectly assumed. The iterations stop 
    % if the enthalpy flag is correct.
    [H,W,T,phi,S,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag] = enthalpy_unsat_flags_iterations...
        (ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflagold,meltold,unsatedgeold,dd,pp,dt);

    % Keep the solutions calculated up until this point so that we can
    % revert to them if we try the other SEB option (melt flag) but find 
    % that that doesn't work.
    Hmid = H;
    msurfmid = msurf;
    Tsurfmid = Tsurf;
    Wmid = W;
    Tmid = T;
    phimid = phi;
    Smid = S;
    hdotmid = hdot;
    pwmid = pw;
    qmid = q;
    roffmid = roff;
   

    %% Step 2: switch surface BC if necessary (and then redo the timestepping)
    % (ie. if we assumed surface melting and got msurf<=0 or if we 
    % assumed no surface melting and got Tsurf>=0)
    %
    % Note: it is important to allow the SEB option to switch if we
    % assumed melt and found msurf=0 or assumed no melt and found
    % Tsurf=0. Both of these situations could be correct, but indicates
    % that we should test the other option too, to prevent the surface
    % getting 'stuck' in a melt or no melt state when really it wants
    % to switch.
    %
    % Again, the enthalpy flag is iterated until it is correct (or until a
    % max number of iterations is reached).

    Tsurftol = 1e-5;

    if and(meltold==1,msurf<=0)        % If melt was assumed but we calculated a negative or 0 melt rate...
        if debug_print == 1
            fprintf([newline,'Swicth to no melt BC (first switch).',newline])
        end
        melt = 0;
        [H,W,T,phi,S,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag] = enthalpy_unsat_flags_iterations...
        (ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflagold,melt,unsatedgeold,dd,pp,dt);    % ... redo the timestep assuming no melt.
        
        if Tsurf>0      % If the new assumption of no melt was also incorrect (Tsurf>0)...
            % ... switch back to the original SEB option.
            melt = 1;
            H = Hmid;
            msurf = msurfmid;
            Tsurf = Tsurfmid;
            W = Wmid;
            T = Tmid;
            phi = phimid;
            S = Smid;
            hdot = hdotmid;
            pw = pwmid;
            q = qmid;
            roff = roffmid;
            if debug_print == 1
                fprintf([newline,'Revert to melt BC.',newline])
            end
            if msurf<0      % If both melt flags give unphysical solutions, keep the second attempt as the solution and set msurf=0.
                if debug_print == 1
                    fprintf([newline,'PROBLEM! Both BCs are unphysical. msurf = ',num2str(msurf),newline]);
                end
                msurf = 0;
            end
        end
            
    elseif and(meltold==0,Tsurf>=0)    % If no melt was assumed but we calculated a positive or 0 surface temp...
        if debug_print == 1
            fprintf([newline,'Swicth to melt BC (first switch).',newline])
        end
            melt = 1;
        [H,W,T,phi,S,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag] = enthalpy_unsat_flags_iterations...
        (ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflagold,melt,unsatedgeold,dd,pp,dt);    % ... redo the timestep assuming melt.

        if msurf<0      % If the new assumption of melt was also incorrect (msurf<0)...
            % ... switch back to the original SEB option.
            melt = 0;
            H = Hmid;
            msurf = msurfmid;
            Tsurf = Tsurfmid;
            W = Wmid;
            T = Tmid;
            phi = phimid;
            S = Smid;
            hdot = hdotmid;
            pw = pwmid;
            q = qmid;
            roff = roffmid;
            if debug_print == 1
                fprintf([newline,'Revert to no melt BC.',newline])
            end
            if Tsurf>0  % If both melt flags give unphysical solutions, keep the second attempt as the solution and set Tsurf=0.
                if debug_print == 1
                    fprintf([newline,'PROBLEM! Both BCs are unphysical. Tsurf = ',num2str(Tsurf),newline]);
                end
                Tsurf = 0;
            end
        end
    end
