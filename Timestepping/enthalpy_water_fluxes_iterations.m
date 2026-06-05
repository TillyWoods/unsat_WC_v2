function [H,W,T,phi,S,Tsurf,msurf,hdot,pw,q,roff,unsatedge,Hflag] = all_iterations_HFWq_new_iter_no_F_it(ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflag,melt,unsatedge,dd,pp,dt)
    % 26-08-2024 - modified to calculate F once at the start of each timestep,
    % rather than inside iterations

    debug_print = 0;    % Set to 1 to print debugging output

    W = ss_old.W;
    phi = ss_old.phi;
    S = ss_old.S;
    % F = ss_old.F;
    q = ss_old.q;
    % a = ss_old.a;
    H = ss_old.H;
    T = ss_old.T;
    n_iterations = 3; %3; %10;  % MUST BE ODD! Otherwise there are jumps in q. (6/6/23)
    for i = 1:n_iterations
        % fprintf([newline,'Iteration ',num2str(i),' out of ',num2str(n_iterations),newline]);
        H1 = H;
        W1 = W;
        phi1 = phi;
        T1 = T;
        S1 = S;
        q1 = q;
        unsatedge1 = unsatedge;
        Hflag1 = Hflag;
        %% Calculate H, W, F from the solution from the previous iteration

        % unsatedge = 0*ones(dd.N+1,1);

        % [H,W,T,phi,S,F,Fplus,Fminus,a,Tsurf,msurf,hdot] = one_iteration_HFW...
            % (ss_old,W,phi,S,F,q,a,alpha,r,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflag,melt,unsatedge,dd,pp,dt);
        % [H,W,T,phi,S,F,Fplus,Fminus,a,Tsurf,msurf,hdot] = one_iteration_HFW_ANALYTICAL_F...
        %     (ss_old,W,phi,S,F,q,a,alpha,r,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflag,melt,unsatedge,dd,pp,dt);
        [H,W,T,phi,S,Tsurf,msurf,hdot] = one_iteration_HW_new_iter_no_F_it...
            (ss_old,a,F,Rs,Rr,Qsi,Q0,nuhat,Tinf,Hflag,melt,unsatedge,dd,pp,dt);

        %% Calculate q from new phi,S just calculated
        % The new flag unsatedgenew (which identifies the edges on which the
        % unsaturated form of the flux should be used) is also calculated.
        % [pw,q,roff,~,unsatedgenew] = water_pressure_and_flux_unsat_flag(phi,S,W,hdot,pp,dd,Rs,Rr);
        [pw,q,roff,~,unsatedgenew] = water_pressure_and_flux_simpler_pw0(phi,S,W,hdot,pp,dd,Rs,Rr);

        %% Update the unsat flag,
        % to be whatever is calculated from the solution of this iteration
        unsatedge_diff = sum(abs(unsatedge-unsatedgenew));
        unsatedge = unsatedgenew;
    
        %% Update the enthalpy flag,
        % to be whatever it is suggested by the solution at this iteration
        [Hflag_new,Hflag_diff] = update_H_flag(H,W,pp,Hflag);
        Hflag = Hflag_new;

        %% Plotting for debugging
        % figure(1)
        % clf;
        % width = 1;
        % hold on
        % p1 = plot(phi,dd.Z,'k-','linewidth',width);
        % p2 = plot(S,dd.Z,'b-','linewidth',width);
        % p3 = plot(10*q,dd.Zedge,'g-','linewidth',width);
        % p4 = plot(pw,dd.Z,'r-','linewidth',width);
        % p5 = plot(T,dd.Z,'m-','linewidth',width);
        % p6 = plot(F,dd.Zedge,'c-','linewidth',width);
        % % p7 = plot([0,1],[Zm,Zm],'k--','linewidth',width);
        % p8 = plot(H,dd.Z,'linewidth',2*width);
        % p9 = plot(W,dd.Z,'linewidth',2*width);
        % set(gca,'YDir','reverse')
        % legend([p1,p2,p3,p4,p5,p6,p8,p9],{'$\phi$','$S$','$10q$','$p_{w}$','$T$','$F$','$H$','$W$'},'interpreter','latex','location','southeast');
        % title(['Iteration ',num2str(i),' out of ',num2str(n_iterations)],'interpreter','latex');
        % ylabel('$Z$','interpreter','latex');
        % xlim([-0.2,1.2]) %xlim([-0.5,0.1]);
        % ylim([0,2]) %ylim([0,3])
        % hold off
        % shg
        
        %% Stop the iterations if both the assumed enthalpy flag and the assumed unsat flag were correct AND the past 2 iterations are within a tolerance
        dif = norm(H1-H) + norm(W1-W); %norm(phi1-phi)+norm(S1-S)+norm(T1-T)+norm(q1-q); % difference from previous iterate
        tol = 1e-3; %1e-6;
        if debug_print == 1
            if and(Hflag_diff==0,unsatedge_diff==0) %and(and(Hflag_diff==0,unsatedge_diff==0),dif<tol)
                fprintf([newline,'Enthalpy flag AND unsat flag correct after ',num2str(i),' out of ',num2str(n_iterations),' iterations. Iterations stopped.',newline])
                break
            else
                fprintf([newline,'Iteration ',num2str(i),' out of ',num2str(n_iterations),': Hflag_diff=',num2str(Hflag_diff),', unsatedge_diff=',num2str(unsatedge_diff),newline]);
                if and(i==n_iterations,Hflag_diff>1)
                    fprintf(['Flags not converged after all iterations, and Hflag_diff=',num2str(Hflag_diff),'!'])
                end
                % if i>7
                %     i
                % end
            end
        end

        % if dif<tol
        %     fprintf([newline,'Iterates converge after ',num2str(i),' out of ',num2str(n_iterations),' iterations. Iterations stopped.',newline])
        %     break
        % else
        %     fprintf([newline,'Iteration ',num2str(i),' out of ',num2str(n_iterations),': Hflag_diff=',num2str(Hflag_diff),', unsatedge_diff=',num2str(unsatedge_diff),', dif=',num2str(dif,3),newline]);
        % end
        if or(sum(W>1)>1,sum(q>1)>1)
            q(3);
        end

    end
end
