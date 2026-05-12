function flag = H_flag(H,W,pp)
    % 13th Dec 2023
    % Calculates the flag that indicates which of the 3 enthlapy cases we
    % are in:
    % flag = 1 if no water
    % flag = 2 if both ice and water
    % flag = 3 if no ice
    % (The value of flag increases with increasing temperature.)

    no_water = H<0;
    no_ice = H>pp.Stefan*W;
    both = and(~no_water,~no_ice);

    flag = 1*no_water + 2*both + 3*no_ice;
end