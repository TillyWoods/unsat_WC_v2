function [H,W,T,phi,S,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag] = enthalpy_unsat_flags_iterations(ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflag,melt,unsatedge,dd,pp,dt)
    % Iterate the enthalpy and water content solve with the water fluxes
    % solve, updating the enthlapy and unsat flags at each interation
    % (which determine in which regions the conductive heat fluxes exist
    % and on which edges the water flux should be included implicity in the
    % enthlpay solver) until both flags agree with the previous iteration,
    % or the maximum number of iterations is reached.

    debug_print = 0;    % Set to 1 to print debugging output

    % Start with the solution from the previous timestep
    W = ss_old.W;           % total water content
    phi = ss_old.phi;       % porosity
    S = ss_old.S;           % saturation
    q = ss_old.q;           % water flux
    H = ss_old.H;           % enthalpy
    T = ss_old.T;           % temperature

    
    % Maximum number of iteration
    max_iterations = 3;  % MUST BE ODD! Otherwise there are jumps in q.
    
    for i = 1:max_iterations
        %% Calculate H, W, F from the solution from the previous iteration
        % The enthalpy flag is also calculated (based on the enthalpy flag
        % from the previous iteration and the new enthalpy solution). This
        % determines where the heat condution term in the enthalpy
        % equation is included (i.e. in regions in which the enthalpy will
        % correspond to temperature).
        [H,W,T,phi,S,Tsurf,msurf,hdot] = enthalpy_water_content_solve...
            (ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflag,melt,unsatedge,dd,pp,dt);

        %% Calculate q from new phi,S just calculated
        % The new flag unsatedgenew (which = 1 on deges that are between 2 
        % unsaturated cells) is also calculated. This determines the cell
        % edges on which the flux q is included implicitly in the enthalpy
        % and water content solve in the next iteration.
        [pw,q,roff,~,unsatedgenew] = water_fluxes_solve(phi,S,W,hdot,pp,dd,Rs,Rr);

        %% Update the unsat flag,
        % to be whatever is calculated from the solution of this iteration
        unsatedge_diff = sum(abs(unsatedge-unsatedgenew));
        unsatedge = unsatedgenew;
    
        %% Update the enthalpy flag,
        % to be whatever it is suggested by the solution at this iteration
        [Hflag_new,Hflag_diff] = update_H_flag(H,W,pp,Hflag);
        Hflag = Hflag_new;
        
        %% Stop the iterations if both the assumed enthalpy flag and the assumed unsat flag were correct
        if and(Hflag_diff==0,unsatedge_diff==0)
            if debug_print == 1
                fprintf([newline,'Enthalpy flag AND unsat flag correct after ',num2str(i),' out of ',num2str(max_iterations),' iterations. Iterations stopped.',newline])
            end
            break
        else
            if debug_print == 1
                fprintf([newline,'Iteration ',num2str(i),' out of ',num2str(max_iterations),': Hflag_diff=',num2str(Hflag_diff),', unsatedge_diff=',num2str(unsatedge_diff),newline]);
                if and(i==max_iterations,Hflag_diff>1)
                    fprintf(['Flags not converged after all iterations, and Hflag_diff=',num2str(Hflag_diff),'!'])
                end
            end
        end

    end
end