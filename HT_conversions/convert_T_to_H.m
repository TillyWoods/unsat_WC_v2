function [H,W] = convert_T_to_H(T,phi,S,pp)
    % 27th October 2022
    W = 1 - phi + phi.*S;               % total water content
    H = W.*T + pp.Stefan*phi.*S;        % enthalpy
end