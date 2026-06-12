function Zedge = nonuniform_Zedge_3_parts(Zend,Zvfine,dZedgevF,Zfine,dZedgeF,dZedgeC,NT)
    % Creates a non-uniform mesh of the cell edge positions, with a very 
    % fine region near the surface, a fine region below, and a coarse 
    % region at depth, and linear transitions between the regions.
    %
    % Zend = end of domain (which starts at Z=0)
    % Zvfine = end of very fine mesh region
    % dZedgevF = very fine cell width
    % Zfine = end of fine mesh region
    % dZedgeF = fine cell width
    % dZedgeC = coarse cell width
    % NT = number of grid cells in transition regions from very fine to
    % fine and from fine to coarse

    NvF = ceil(Zvfine/dZedgevF);                                                    % number of very fine grid cells
    ZedgevF = linspace(0,NvF*dZedgevF,NvF+1)';                                          % very fine mesh
    dZedgeT1 = dZedgevF + (dZedgeF-dZedgevF)*(linspace(NvF+1,NvF+NT,NT)'-NvF)/(NT+1);     % cell widths in transition region
    ZedgeT1 = zeros(NT,1);                                                           % transition mesh from very fine to fine
    for i = 1:NT
        ZedgeT1(i) = ZedgevF(end) + sum(dZedgeT1(1:i));
    end
    
    NF = ceil((Zfine-ZedgeT1(end))/dZedgeF);                                                       % number of fine grid cells
    ZedgeF = ZedgeT1(end)+dZedgeF + linspace(0,(NF-1)*dZedgeF,NF)';                                          % fine mesh
    dZedgeT2 = dZedgeF + (dZedgeC-dZedgeF)*(linspace(NF+1,NF+NT,NT)'-NF)/(NT+1);     % cell widths in transition region
    ZedgeT2 = zeros(NT,1);                                                           % transition mesh fine to coarse
    for i = 1:NT
        ZedgeT2(i) = ZedgeF(end) + sum(dZedgeT2(1:i));
    end

    NC = ceil((Zend-ZedgeT2(end))/dZedgeC);                                          % number of coarse grid cells
    dZedgeCnew = (Zend-ZedgeT2(end))/NC;                                             % adjusted coarse cell widths so there is an integer number of cells
    ZedgeC = linspace(ZedgeT2(end)+dZedgeCnew,Zend,NC)';                            % coarse mesh
    
    Zedge = [ZedgevF;ZedgeT1;ZedgeF;ZedgeT2;ZedgeC];                                                 % full mesh
end
