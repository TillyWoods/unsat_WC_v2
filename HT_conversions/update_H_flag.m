function [new_flag,wrong] = update_H_flag(H,W,pp,flag)
    % 13th Dec 2023
    % Updates the values of the enthalpy flag based on the solution H,W and
    % the flag we assumed for calculating that solution. If the assumed
    % flag is inconsistent with the solution, it is updated (see below).
    %
    % NOTE: we do not simply update the flag to be that calculated from the
    % solution! We have to be more careful to ensure that a grid cell can't
    % jump straight from only-ice to only-water or vice versa. It must go
    % through a both-ice-and-water state first.
    %
    % - wrong = the number of entries of the flag that were incorrectly
    % assumed.
    
    % What values of the flag mean:
    % - flag = 1 for ice-only (T<0, phi=0, H<0)
    % - flag = 2 for both ice and water (T=0, 0<phi<1, 0<=H<=Stefan*W)
    % - flag = 3 for water-only (T>0, phi=1, H>Stefan*W)
    
    new_flag = flag;
    sol_flag = H_flag(H,W,pp);  % flag calculated using the solution, to tell us what state the solution is in

    % If the state that the solution is in (only ice, both ice and water,
    % or only water) is inconsistent with what we assumed (flag), update 
    % the value of the flag as follows:
    error_1to2 = and(flag==1,~(sol_flag==1));     % If assumed only ice but found T>0 (H>0, sol_flag=2,3)...
    new_flag(error_1to2) = 2;                   % ... try assuming there's both ice and water
    error_2to1 = and(flag==2,sol_flag==1);        % If assumed both ice and water but found phi<0 (H<0, sol_flag=1)...
    new_flag(error_2to1) = 1;                   % ... try assuming there's only ice
    error_2to3 = and(flag==2,sol_flag==3);        % If assumed both ice and water but found phi>1 (H>Stefan*W, sol_flag=3)...
    new_flag(error_2to3) = 3;                   % ... try assuming there's only water
    error_3to2 = and(flag==3,~(sol_flag==3));     % If assumed only water but found T<0 (H<0, sol_flag=1)...
    new_flag(error_3to2) = 2;                   % ... try assuming there's both ice and water

    % Number of entries which where inconsistent (i.e. assuming flag didn't
    % agree with solution):
    wrong = sum(error_1to2) + sum(error_2to1) + sum(error_2to3) + sum(error_3to2);
end