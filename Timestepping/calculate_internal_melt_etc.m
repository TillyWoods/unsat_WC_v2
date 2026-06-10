function ss = calculate_internal_melt_etc(ss,ss_old,pp,dd,dt)
    % Calculate additional useful quantities at the end of the timestep:
    % 
    % - mint - subsurface melt rate
    % - Zm - bottom of weathering crust (measured as distance below the
    % surface)
    % - Zf - top of weathering crust/freezing from (measured as distance 
    % below the surface)
    % - crust - logical that is 1 where weathering crust is present
    % - h - surface height
    % - qd - lateral drainage
    % - musurf_cum - sumulative surface melt rate
    % - mint_cum - cumulative subsurface melt rate NOT INCLUDING FREEZING
    % - freeze_cum - cumulative subsurface FREEZING
    % - roff_cum - cumulative surface runoff
    % - qd_cum - cumulative lateral drainage

    
    %% Internal melt rate, mint
    ss.mint = (ss.phi-ss_old.phi)/dt ...
        + dd.Dflux*(dd.cell_below*min(ss.hdot,0)...
        +dd.cell_above*max(ss.hdot,0))*ss.phi;

    %% Locate bottom of the weathering crust, Zm...
    % ...using zero contour matrix of enthalpy. (Two ss.H are used to form 
    % a matrix so contourc works).
    M = contourc([1,2],dd.Z,[ss.H,ss.H],[0,0]);
    % Zm is calculated as the max value of the contour. The contour values 
    % are every 3rd entry f the 2nd row, starting from the second entry.
    ss.Zm = max(M(2,2:3:end));                      
    ss.Zm(isempty(ss.Zm)) = 0;

    %% Locate top of the weathering crust (/freezing front), Zf 
    % Zf (freezing front/top of weathering crust) is calculated as the min 
    % value of the contour. The contour values are every 3rd entry f the 
    % 2nd row, starting from the second entry.
    ss.Zf = min(M(2,2:3:end));                      
    ss.Zf(isempty(ss.Zf)) = 0;      % if there is no reeezing front...
    ss.Zf(ss.Zf==ss.Zm) = 0;        % ... or if it Zm is found as the only contour, then set Zf to be at the surface, i.e. 0 in the Z frame

    ss.crust = ss.Zm - ss.Zf > 1e-12;      % logical that is 1 if there is a weathering crust and 0 is not

    %% Surface height, h
    ss.h = ss_old.h + dt*ss.hdot;

    %% Lateral drainage, qd
    ss.qd = pp.K*pp.Kd*pp.perm(ss.phi).*pp.krw(ss.S);

    %% Cumulative surface and subsurface melting and runoff, 
    % and subsurface freezing
    ss.msurf_cum = ss_old.msurf_cum + dt*ss.msurf;          % cumulative surface melt
    ss.mint_cum = ss_old.mint_cum + dt*max(ss.mint,0);      % cumulative subsurface MELT, NOT INCLUDING FREEZING
    ss.freeze_cum = ss_old.freeze_cum + dt*min(ss.mint,0);  % cumulative subsurface FREEZING
    ss.roff_cum = ss_old.roff_cum + dt*ss.roff;             % cumulative surface runoff
    ss.qd_cum = ss_old.qd_cum + dt*ss.qd;                   % cumulative lateral drainage
end