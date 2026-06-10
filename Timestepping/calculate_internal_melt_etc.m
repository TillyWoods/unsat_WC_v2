function ss = calculate_internal_melt_etc(ss,ss_old,pp,dd,dt)
    % Calculate additional useful quantities at the end of the timestep
    % (internal melt rate etc.)

    ss.mint = (ss.phi-ss_old.phi)/dt + dd.Dflux*(dd.cell_below*min(ss.hdot,0)+dd.cell_above*max(ss.hdot,0))*ss.phi;    % internal melt rate     % EXTRA hdot (INCORRECT!) REMOVED 08/04/2024
    % ss.Zm = dd.Z(find(ss.phi>0,1,'last'));          % last entry of Z that is non-zero (i.e. phi=0 everywhere below)
    M = contourc([1,2],dd.Z,[ss.H,ss.H],[0,0]);     % zero contour matrix of enthalpy. Two ss.H are used to form a matrix so contourc works.
    ss.Zm = max(M(2,2:3:end));                      % Zm is calculated as the max value of the contour. The contour values are every 3rd entry f the 2nd row, starting from the second entry.
    ss.Zm(isempty(ss.Zm)) = 0;
    ss.h = ss_old.h + dt*ss.hdot;
    ss.Zf = min(M(2,2:3:end));                      % Zf (freezing front/top fof weathering crust) is calculated as the min value of the contour. The contour values are every 3rd entry f the 2nd row, starting from the second entry.
    % if isempty(ss.Zf)
    %     ss.Zf
    % end
    
    ss.Zf(isempty(ss.Zf)) = 0;      % if there is no reeezing front...
    % ss.Zf(isempty(ss.Zf)) = 0;
    ss.Zf(ss.Zf==ss.Zm) = 0;        % ... or if it Zm is found as the only contour, then set Zf to be at the surface, i.e. 0 in the Z frame

    ss.crust = ss.Zm - ss.Zf > 1e-12;      % logical that is 1 if there is a weathering crust and 0 is not

    ss.msurf_cum = ss_old.msurf_cum + dt*ss.msurf;  % cumulative surface melt (ADDED 29/01/2024)
    ss.mint_cum = ss_old.mint_cum + dt*max(ss.mint,0);     % cumulative internal MELT, NOT INCLUDING FREEZING (ADDED 29/01/2024)
    ss.freeze_cum = ss_old.freeze_cum + dt*min(ss.mint,0);  % cumulative internal FREEZING (ADDED 09/08//2024)
    ss.roff_cum = ss_old.roff_cum + dt*ss.roff;  % cumulative surface runoff (ADDED 17/09/2024)

    ss.qd = pp.K*pp.Kd*pp.perm(ss.phi).*pp.krw(ss.S);   % lateral drainage (ADDED 13/03/2024, with lengthscale over which drainage occurs absorbed into the dimensionless form of qd)
    ss.qd_cum = ss_old.qd_cum + dt*ss.qd;                   % cumulative lateral drainage (ADDED 13/03/2024 - CORRECTED 18/09/2024)
end