%% Periodic switching experiments - 30/05/2025
%
% A periodic version of the switching experiment in the unsaturated
% chapter of my thesis.
%
% Using the implicit, unsaturated code (with lateral drainage), but for
% runs both with and without lateral drainage.
%
% The forcing goes periodically between 3 sections of constant forcing
% values, detailed below.

%%
addpath(genpath('.'))   % Add all subfolders to path

Kd_vals = 50; %[0,3];

% Use no drainage solution after 10 periods as initial condition, if the
% logical below is set to 1 (added 14/08/2025). Otherwise, use steadily
% melting initial condition
%load("wc_periodic_switch_Kd_0_14-Aug-2025_15_59_12.mat")
% load("wc_periodic_switch_Kd_3_19-Aug-2025_13_04_20.mat")   % 3 drainage periods run after the 10 no-drainage periods above (19/08/2025)
% load("wc_periodic_switch_Kd_0_19-Aug-2025_11_33_55.mat")   % 3 no drainage periods run after the 10 no-drainage periods above (19/08/2025)
% ss_init = ss_save(end);
use_10period_init = 0;


for i = 1:length(Kd_vals)
%%% Saturated or unsaturated?
Kd = Kd_vals(i);
% Use smaller timesteps if water flow is allowed
if Kd == 0
    dt = 0.25*1e-3; %1e-2; %0.25*1e-3;
    gap = 4*1e2; %10;
else
    dt = 0.25*1e-3;
    gap = 4*1e2;
end

%%% Discretisation
dd = discretise_isothermal();
dd.plot = 0;
dd.gap = gap; %1e2; %10;
dd.dt = dt;
dd.steady_stop = 0;
dd.steady_tol = 1e-4; %1e-4;

%%% Parameters and forcings
dim = dimensional_parameters();
% Values used for steadily melting initial condition
Qsi = 100; %200; %0;
Q0 = 50; %-20; %20; %100; 

% First section:
Qsi1 = 100;
Q01 = 50;
gap1 = 20;      % duration of section

% Second section:
Qsi2 = 0;
Q02 = -50;
gap2 = 20;      % duration of section

% Third section:
Qsi3 = 0;
Q03 = 100;
gap3 = 20;      % duration of section

period = gap1 + gap2 + gap3;

dim.Qsi = @(t) Qsi1*(mod(t,period*dim.td)<gap1*dim.td)...
    +Qsi2*and(mod(t,period*dim.td)>=gap1*dim.td,mod(t,period*dim.td)<=(gap1+gap2)*dim.td)...
    +Qsi3*(mod(t,period*dim.td)>(gap1+gap2)*dim.td);
dim.Q0 = @(t) Q01*(mod(t,period*dim.td)<gap1*dim.td)...
    +Q02*and(mod(t,period*dim.td)>=gap1*dim.td,mod(t,period*dim.td)<=(gap1+gap2)*dim.td)...
    +Q03*(mod(t,period*dim.td)>(gap1+gap2)*dim.td);

scale = scalings(dim);
pp = dimensionless_parameters(dim,scale);

% Drainage or not?
pp.Kd = Kd;

% pp.B = 1e4;

pp.tolphi0 = 1e-4;

%%% Initial condition
if use_10period_init==1
    % Use end of 10 period no-drainage run as initial condition (added
    % 14/08/2025)
    ss0 = ss_init;
else
    % Use steadily melting as initial condition
    ss0 = steadily_melting_initial_cond(Qsi/scale.F,Q0/scale.F,pp,dd);  % Dimensionless Qsi and Q0 in here
    % ss0.phi = 0*dd.Z.^0;
    ss0.S = 1*dd.Z.^0; %0.5*dd.Z.^0;        % fully saturated initially
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp);
end

tspan =  [0,5]; %20*period]; %[0,10*period]; %[0,1]; %[0,100]; %[0,0.1];


%% Timestepping
tic
ss = timestepping_unsat_implicit(ss0,pp,dd,tspan);
toc

%% Save data
save_time = 240; %inf;      % save the last save_time amount of time of data
save_ind = ss(end).t-[ss.t]<save_time;
% save_ind = 1:40:length(ss);
ss_save = ss(save_ind);
tsave = [ss_save.t];

save(["wc_periodic_switch_Kd_"+num2str(pp.Kd,2)+"_"+strrep(strrep(string(datetime),' ','_'),':','_')+".mat"],"ss_save","pp","dd","scale","dim")

end
%% Plotting
h = [ss.h];
phi = [ss.phi];
S = [ss.S];
T = [ss.T];
t = [ss.t];
hdot = [ss.hdot];
roff = [ss.roff];
msurf = [ss.msurf];
Zm = [ss.Zm];
H = [ss.H];
Zf = [ss.Zf];

[hgrid,Zgrid] = meshgrid(h,dd.Z);
[tgrid,Zgrid] = meshgrid(t,dd.Z);
zgrid = hgrid - Zgrid;

[hedge,Zedge] = meshgrid(h,dd.Zedge);
[tedge,Zedge] = meshgrid(t,dd.Zedge);
zedge = hedge - Zedge;

ztop = 0;
zlimit = -3;
Ztop = 0;
Zlimit = 0.5;
tstart = 0; %tspan(1);
tend = tspan(2);
fixed_surface = 0;      % 1 %%

width = 1;

depth_grid = fixed_surface*Zgrid + (1-fixed_surface)*zgrid;

TphiS_plot = figure(3);
box on
set(gcf,'units','centimeters','paperpositionmode','auto','position',[2 2 9 10.6]);
tiledlayout(4,1,'TileSpacing','compact')

colormap(parula(256));

ax(1) = nexttile;
hold on
surf(tgrid,depth_grid,T*scale.T,'LineStyle','none')
view(0,90);
% contour(tgrid,depth_grid,T*scale.T,[-1e-12,-1e-12],'color','k','linewidth',1.5,'linestyle','--')
% contour(tgrid,depth_grid,H*scale.H,[0,0],'color','k','linewidth',1.5,'linestyle','-')
xlim([min(t),max(t)])
set(gca,'XTick',[])
c1 = colorbar;
% caxis([-0.1,0])
ylabel(c1,'$T$ ($^{\mathrm{o}}$C)','interpreter','latex')
if fixed_surface == 1
    plot3(t(Zm>0),Zm(Zm>0),max(T,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    plot3(t(Zf>0),Zf(Zf>0),max(T,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel('$Z$ (m)','interpreter','latex')
    set(gca,'YDir','reverse')
    caxis([-(-min(scale.T*T(Zgrid<Zlimit),[],'all')),0])
else
    plot3(t(Zm>0),h(Zm>0)-Zm(Zm>0),max(T,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    plot3(t(Zf>0),h(Zf>0)-Zf(Zf>0),max(T,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel('$z$ (m)','interpreter','latex')
    plot3(t,h,max(T,[],'all')*t.^0,'k-','linewidth',width)
    caxis([-(-min(scale.T*T(zgrid>zlimit),[],'all')),0])
end
% caxis([-10,0])

ax(2) = nexttile;
hold on
surf(tgrid,depth_grid,phi,'LineStyle','none')
view(0,90);
% contour(tgrid,depth_grid,phi,[1e-10,1e-10],'color','k','linewidth',1.5,'linestyle','--')
% contour(tgrid,depth_grid,H*scale.H,[0,0],'color','k','linewidth',1.5,'linestyle','-')
xlim([min(t),max(t)])
% title('$\phi$','interpreter','latex')
set(gca,'XTick',[])
c2 = colorbar;
% caxis([0,max(max(phi(and(t>=tstart,t<=tend))))])
caxis([0,0.2])
ylabel(c2,'$\phi$','interpreter','latex')
if fixed_surface == 1
    plot3(t(Zm>0),Zm(Zm>0),max(phi,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    plot3(t(Zf>0),Zf(Zf>0),max(phi,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel('$Z$ (m)','interpreter','latex')
    set(gca,'YDir','reverse')
else
    plot3(t(Zm>0),h(Zm>0)-Zm(Zm>0),max(phi,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    plot3(t(Zf>0),h(Zf>0)-Zf(Zf>0),max(phi,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel('$z$ (m)','interpreter','latex')
    plot3(t,h,max(phi,[],'all')*t.^0,'k-','linewidth',width)
end

% ax(3) = nexttile;
% surf(tgrid,zgrid,S,'LineStyle','none')
% view(0,90);
% xlim([min(t),max(t)])
% xlabel('$t$ (day)','interpreter','latex')
% ylabel('$z$ (m)','interpreter','latex')
% title('$S$','interpreter','latex')
% colorbar
% caxis([0,1])

% linkaxes(ax)

% ylim([-inf,0])


ax(3) = nexttile;
hold on
plot(t,dim.Qsi(t*scale.t),'r-','displayname','$Q_\mathrm{si}$')
plot(t,dim.Q0(t*scale.t),'b-','displayname','$Q_0$')
set(gca,'XTick',[])
% xlabel('$t$ (day)','interpreter','latex')
ylabel('(W m$^{-2}$)')
legend('show','interpreter','latex','position',[0.77,0.36,0.2,0.1]) %,[0.75,0.2,0.2,0.1])

ax(4) = nexttile;
hold on
yyaxis left
plot(t,-hdot,'displayname','$-\dot{h}$')
% plot(t,roff,'r--','displayname','$r_\mathrm{off}$')
plot(t,msurf,'displayname','$m_\mathrm{surf}$','linestyle',':')
ylabel('(m d$^{-1}$)')
% ylim([5.6e-3,6.4e-3])
yyaxis right
plot(t,Zm,'displayname','$Z_\mathrm{m}$')
plot(t,mean(Zm(and(t>=tstart,t<=tend)))*Zm.^0,'displayname','Mean','linestyle','--')
% plot(t,Zmsteady*Zm.^0,'displayname','Steady','linestyle',':')
ylabel('(m)')
% ylim([0,0.4])
% xlabel('$t$ (day)','interpreter','latex')
leg1 = legend('show','interpreter','latex','position',[0.8,0.13-0.12,0.23-0.05,0.14]);
leg1.ItemTokenSize = [15,1];

linkaxes(ax,'x')
linkaxes([ax(1),ax(2)],'y')
if fixed_surface == 1
    ylim(ax(1),[Ztop,Zlimit])
else
    ylim(ax(1),[zlimit,ztop])
end
xlim([tstart,tend])
