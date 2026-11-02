%% Figure 8: Melt and runoff components
% 08/10/2025
% Plotting the surface and subsurface melt and runoff (and freezing)
% components over 1 period of the periodic switching run.

set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1.5)

addpath(genpath('./..'))       % Add all subfolders to path

% Directory to save plots into
output_dir = ".";

plot_extras = 0;

% Colours
cmap8 = colblind8();
mint_col = cmap8(5,:);
msurf_col = cmap8(3,:);
qd_col = cmap8(1,:);
freeze_col = cmap8(6,:);
roff_col = cmap8(4,:);
loss_col = 0.8*[1 1 1];

Z = ss(1).Z;            % Z coordinates
Zedge = ss(1).Zedge;
t = [ss.t] - (ss(end).t - 62);  
msurf = [ss.msurf];     % surface melt rate
mint = [ss.mint];       % subsurface melt rate (at each Z)
roff = [ss.roff];       % surface runoff rate
qd = [ss.qd];           % subsurface lateral drainage flux (at each Z)
meltin = max(mint,0);   % subsurface melt rate (EXCLUDING freezing) (at each Z)
freeze = -min(mint,0);   % subsurface freezing rate (at each Z)
q = [ss.q];

top_id = 2;

qdZ = sum(qd(top_id:end,:).*diff(Zedge(top_id:end))); % trapz(Z,qd);          % total subsurface lateral drainage rate (integrated over Z)
meltinZ = sum(meltin(top_id:end,:).*diff(Zedge(top_id:end))); %trapz(Z,meltin);  % total subsurface melt rate (integrated over Z)
freezeZ = sum(freeze(top_id:end,:).*diff(Zedge(top_id:end))); %trapz(Z,freeze);  % total subsurface freezing rate (integrated over Z)
mintZ = trapz(Z,mint);

% Recalculate all the internal melt rates etc. manually for debugging
% (02/06/2026)
mint_manual = ([ss.phi] - [ss(1).phi,[ss(1:end-1).phi]])/dd.dt/10; % [zeros(1,length([ss.t]));diff([ss.phi],1)/dd.dt];       % mint calculated as (phi-phi_old)/dt
melt_manual = max(mint_manual,0);
freeze_manual = -min(mint_manual,0);
melt_manualZ = trapz(Z,melt_manual);
freeze_manualZ = trapz(Z,freeze_manual);
mint_manualZ = trapz(Z,mint_manual);

% Where to start and end the period of interest
shift = 2;
tstart = t(end)-60-shift;
tend = t(end) - shift;
istart = find(t>=tstart,1,'first');
iend = find(t>=tend,1,'first');


W = [ss.W];         % total water content
S = [ss.S];
phi = [ss.phi];
q = [ss.q];
hdot = [ss.hdot];
tolS1 = 1e-4; %pp.tolS1;
satt = ((3/2)*W(1,:)-(1/2)*W(2,:) >= 1-tolS1);      % is the top cell (representing surface) saturated?
roffnew = zeros(size(roff));
qsurfsat = 2*q(2,:)-q(3,:);     % linearly extrapolate surface q from the next two
roffsat = -qsurfsat - [hdot(1),hdot(1:end-1)].*W(1,:);    % calculate runoff in saturated cells
roffnew(satt) = roffsat(satt);

% Integrate over time (cumulative)
% Note: the internal melt and freeze, lateral drainage and surface melt can
% be taken from the saved solution (should be correct). But they should be
% shifted to be zero at the start of the period. Cumulative surface
% runoff should be recalculated.
% NOTE: the below might not be completely correct because they are
% calculated using the saved timesteps, which is not all timesteps.#

freeze_cum = [ss.freeze_cum];
meltin_cum = [ss.mint_cum];
qd_cum = [ss.qd_cum];

msurf_cum = [ss.msurf_cum]; %cumtrapz(msurf); %[ss.msurf_cum];
freezeZ_cum = -sum(freeze_cum(top_id:end,:).*diff(Zedge(top_id:end))); %-trapz(Z,[ss.freeze_cum]); %cumtrapz(t,freezeZ);
meltinZ_cum = sum(meltin_cum(top_id:end,:).*diff(Zedge(top_id:end))); % trapz(Z,[ss.mint_cum]); % cumtrapz(t,meltinZ);
qdZ_cum = sum(qd_cum(top_id:end,:).*diff(Zedge(top_id:end))); %trapz(Z,[ss.qd_cum]); % cumtrapz(t,qdZ);
% freezeZ_cum = -sum(freeze_cum.*diff(Zedge)); %-trapz(Z,[ss.freeze_cum]); %cumtrapz(t,freezeZ);
% meltinZ_cum = sum(meltin_cum.*diff(Zedge)); % trapz(Z,[ss.mint_cum]); % cumtrapz(t,meltinZ);
% qdZ_cum = sum(qd_cum.*diff(Zedge)); %trapz(Z,[ss.qd_cum]); % cumtrapz(t,qdZ);
roff_cum = [ss.roff_cum]; %cumtrapz(t,roffnew);

% msurf_cum = cumtrapz(msurf); %[ss.msurf_cum];
% freezeZ_cum = cumtrapz(t,freezeZ);
% meltinZ_cum = cumtrapz(t,meltinZ);
% qdZ_cum = cumtrapz(t,qdZ);
% roff_cum = cumtrapz(t,roffnew);

shift_start = @(x) x - x(istart);   % shift to be zero at the start of period of interest

msurf_cum = shift_start(msurf_cum);
freezeZ_cum = shift_start(freezeZ_cum);
meltinZ_cum = shift_start(meltinZ_cum);
qdZ_cum = shift_start(qdZ_cum);
roff_cum = shift_start(roff_cum);

melt_total = msurf_cum + meltinZ_cum;
runoff_total = roff_cum + qdZ_cum;
removal_total = runoff_total + freezeZ_cum;
net = melt_total - removal_total;      % net mass production - should be zero at end of period

%% Melt and runoff components over time (one period)
f7 = figure(8);
set(gcf,'PaperPositionMode','auto','Units','centimeters','position',[4 4 14 4])
clf
box on
hold on
plot(t,msurf,'DisplayName','Surface melt','Color',msurf_col)
plot(t,meltinZ,'DisplayName','Subsurface melt','Color',mint_col)
plot(t,roff,'DisplayName','Surface runoff','Color',roff_col)
plot(t,qdZ,'DisplayName','Subsurface runoff','Color',qd_col)
plot(t,freezeZ,'DisplayName','Freezing','Color',freeze_col)
%plot(t,freeze_manualZ,'DisplayName','Freezing manual','Color',freeze_col,'LineStyle','--')
%plot(t,melt_manualZ,'DisplayName','Subsurface melt manual','color',mint_col,'linestyle','--')
xline([2,22,42],'linestyle','--','HandleVisibility','off')
%plot(t,q(1,:),'DisplayName','Surface q','LineStyle','--')
% plot(t,W(1,:)/50)
% plot(t,roff)
xlabel('$t$')
ylabel('Rate')
if pp.Kd==0
    leg = legend('show','interpreter','latex','fontsize',9', 'position',[0.75,0.3,0.2,0.5]); %,location','east');
    leg.ItemTokenSize = [10 5];
end
xlim([tstart,tend])
% ylim([0,inf]) %([0,0.032])

print(f7,output_dir+"/melt_runoff_Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-'),'-dpng','-r600')
print(f7,output_dir+"/melt_runoff_Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-'),'-depsc','-r600')
print(f7,output_dir+"/melt_runoff_Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-'),'-dsvg','-r600')
%print(f7,output_dir+"/melt_runoff_Kd_"+num2str(pp.Kd,2)+"_tolphi0_"+strrep(num2str(pp.tolphi0),'.','-')+"_tolS1_"+strrep(num2str(pp.tolS1),'.','-'),'-dpng','-r600')

if plot_extras == 1
    %% Surface lowering, rounoff, surface fluxes etc.
    % New bit 17/10/2025: to test issues with using different conditions to
    % decide whether the surface is saturated for calculating surfacewater flux
    % vs. surface runoff
    
    % For calcualting surface water flux:
    sat = and(S>=1-pp.tolS1,phi>0+pp.tolphi0);
    solid = logical(phi<=pp.tolphi0);
    
    % For calculating surface runoff:
    sat_for_roff = (3/2)*W(1,:)-(1/2)*W(2,:) >= 1 - pp.tolS1;
    
    f1 = figure(11);
    clf
    box on
    hold on
    plot(t,-hdot*30,'DisplayName','-hdot')
    plot(t,q(1,:)*30,'displayname','q')
    plot(t,(-hdot-q(1,:))*30,'displayname','-hdot-q')
    plot(t,roff*30,'DisplayName','roff')
    yline(1)
    plot(t,W(1,:),'DisplayName','W1')
    plot(t,W(2,:),'DisplayName','W2')
    plot(t,(3/2)*W(1,:) - (1/2)*W(2,:),'DisplayName','Wsurf','LineStyle','--')
    plot(t,sat(1,:)+solid(1,:),'m:','DisplayName','sat(1)+solid(1)')
    plot(t,sat_for_roff,'g:','DisplayName','Wsurf>=1-tolS1')
    plot(t,sat(1,:),'b:','DisplayName','sat(1)')
    plot(t,solid(1,:),'r:','DisplayName','solid(1)')
    plot(t,(3/2)*S(1,:) - (1/2)*S(2,:),'DisplayName','Ssurf','LineStyle','--')
    plot(t,phi(1,:),'DisplayName','phi(1)')
    plot(t,(3/2)*phi(1,:)-(1/2)*phi(2,:),'DisplayName','phisurf')
    plot(t,(3/2)*phi(1,:)-(1/2)*phi(2,:) <= pp.tolphi0,'DisplayName','phisurf<=tolphi0','LineStyle',':')
    xlabel('$t$')
    ylabel('$-\dot{h}$')
    xlim([tstart,tend])
    legend('show')
    
    %% Cumulative melt and runoff
    f3 = figure(3);
    clf
    hold on
    box on
    plot(t,melt_total,'DisplayName','Total melt')
    plot(t,runoff_total,'DisplayName','Total runoff')
    plot(t,freezeZ,'DisplayName','Total freeze')
    plot(t,net,'DisplayName','Net mass gain')
    plot(t,msurf_cum,'DisplayName','Surface melt','LineStyle','--')
    plot(t,meltinZ_cum,'DisplayName','Subsurface melt','LineStyle','--')
    plot(t,roff_cum,'DisplayName','Surface runoff','LineStyle','--')
    plot(t,qdZ_cum,'DisplayName','Subsurface runoff','LineStyle','--')
    plot(t,q(1,:),'DisplayName','Surface q','LineStyle','--')
    xlabel('$t$')
    ylabel('Cumulative...')
    xlim([tstart,tend])
    legend('show')
end