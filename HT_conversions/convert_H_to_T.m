function [T,phi,S] = convert_H_to_T(H,W,pp)
    % Converts from enthalpy and total water content to temperature,
    % porosity and saturation. It does this for the 3 cases: ice only,
    % water only, and both ice and water (update in Nov/Dec 2023).

    no_water = H<0;
    no_ice = H>pp.Stefan*W;
    both = and(~no_water,~no_ice);
    
    %%%%
    % Update 13/5/24
    % Tolerances below/above which to set phi to 0 or S to 0 or 1.
    % Previously, I have not SET S or phi to be these values, even if I had
    % used the tolerances when 'dividing by 0'.
    tolphi0 = pp.tolphi0;
    tolS0 = pp.tolS0;
    tolS1 = pp.tolS1;

    tol = tolphi0;
    %%%%

    T = no_water.*H./W + no_ice.*(H./W - pp.Stefan);
    phi = no_water.*(1- W) + both.*(1 - W + H/pp.Stefan) + no_ice; %1 - (1-no_ice).*W + both.*H/pp.Stefan;
    S = both.*H./(pp.Stefan*max(phi,tol)) + no_ice.*W;
    %S(phi<tol) = 0; %1; %0;      % NOT SURE IF THIS SHOULD BE 0 OR 1?

    %%%%
    % Update 13/5/24
    %phi(phi<=tolphi0) = 0;
    %S(S<=tolS0) = 0;
    %S(S>=1-tolS1) = 1;
    %%%%
end