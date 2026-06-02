function [H,W,T,phi,S,F,Fplus,Fminus,a,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag,melt] = HFWq_timestep_new_iter_no_F_it(ss_old,alpha,r,Rs,Rr,Qsi,Q0,nuhat,Tinf,dd,pp,dt)
% 26-08-2024 - modified to calculate F once at the start of each timestep,
% rather than inside iterations

debug_print = 0;    % Set to 1 to print debugging output

meltold = ss_old.melt;
melt = meltold;
Hflagold = ss_old.Hflag;
unsatedgeold = ss_old.unsatedge;

%% Solve for F (using phi, S from the previous timestep, if applicable)
% phihm = harmonic_mean(ss_old.phi);
% Shm = harmonic_mean(ss_old.S);
% [Fplusnew,Fminusnew] = F_solve_second_order(Qsi,alpha(phihm,Shm),r(phihm,Shm),pp.chi,dd);
% a = pp.a;

alpha0 = alpha(1,1);
alpha0 = alpha0(1);
r0 = r(1,1);
r0 = r0(1);
kappa = sqrt(alpha0^2+2*alpha0*r0);
a0 = (1-pp.chi)*(kappa-alpha0)/((1-pp.chi)*kappa+(1+pp.chi)*alpha0);
if pp.chi==1
        Fminus = zeros(size(dd.Zedge));
        Fplus = zeros(size(dd.Zedge));
else
    Fminus = (1-pp.chi*(1-a0))*Qsi*exp(-kappa*dd.Zedge);
    Fplus = a0*Qsi*exp(-kappa*dd.Zedge);
end

F = Fminus - Fplus;
a = a0;
% a = (Fplus(1)/Qsi);
% a(Qsi==0) = a; %pp.a;   % if Qsi=0 (in which case F=0), set a to be the default value



 %% Step 1: carry out the timestep...
    % ... assuming that the melt and enthlapy flags are that of the
    % previous timestep.
    % NOTE: within this, the a number (probably 10 - see function below to 
    % be sure) of iterations are carried out in which the enthalpy flag is 
    % updated each time if it was incorrectly assumed. The iterations stop 
    % if the enthalpy flag is correct.
    [H,W,T,phi,S,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag] = all_iterations_HFWq_new_iter_no_F_it...
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
        [H,W,T,phi,S,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag] = all_iterations_HFWq_new_iter_no_F_it...
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
        [H,W,T,phi,S,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag] = all_iterations_HFWq_new_iter_no_F_it...
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