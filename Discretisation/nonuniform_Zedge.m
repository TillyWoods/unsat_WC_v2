function Zedge = nonuniform_Zedge(Zend,Zfine,dZedgeF,dZedgeC,NT)
    % 2nd November 2022
    % Creates a non-uniform mesh of the cell edge positions, with a fine
    % region near the surface, a coarse region at depth, and a linear
    % transition region inbetween.
    %
    % Zend = end of domain (which starts at Z=0)
    % Zfine = end of fine mesh region
    % dZedgeF = fine cell width
    % dZedgeC = coarse cell width
    % NT = number of grid cells in transition region between fine and coarse mesh
    
    NF = ceil(Zfine/dZedgeF);                                                       % number of fine grid cells
    ZedgeF = linspace(0,NF*dZedgeF,NF+1)';                                          % fine mesh
    dZedgeT = dZedgeF + (dZedgeC-dZedgeF)*(linspace(NF+1,NF+NT,NT)'-NF)/(NT+1);     % cell widths in transition region
    ZedgeT = zeros(NT,1);                                                           % transition mesh
    for i = 1:NT
        ZedgeT(i) = ZedgeF(end) + sum(dZedgeT(1:i));
    end
    NC = ceil((Zend-ZedgeT(end))/dZedgeC);                                          % number of coarse grid cells
    dZedgeCnew = (Zend-ZedgeT(end))/NC;                                             % adjusted coarse cell widths so there is an integer number of cells
    ZedgeC = linspace(ZedgeT(end)+dZedgeCnew,Zend,NC)';                             % coarse mesh
    Zedge = [ZedgeF;ZedgeT;ZedgeC];                                                 % full mesh
end
