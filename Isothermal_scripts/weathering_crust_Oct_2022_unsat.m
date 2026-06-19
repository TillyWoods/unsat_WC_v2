% Weathering crust code for when there is unsaturated flow and lateral
% drainage
% Feb 2024 (modified from October 2022)

% Used for both isothermal water transport test case and with real KAN_M
% forcng data

%%
addpath(genpath('/scratch/woods/Dropbox/Oxford Mathematics DPhil/MATLAB/Weathering crust/Returning to WC Oct 2022'))
addpath(genpath('C:\Users\twood\Dropbox\Oxford Mathematics DPhil\MATLAB\Weathering crust/Returning to WC Oct 2022'))

%% Enter the run numner
% 1 for isothermal water transport
% 2 for KAN_M PROMICE data
run_num = 1;

if run_num == 1
    %% Isothermal water transport
    % Discretisation
    dd = discretise_isothermal();
    dd.plot = 0; %1;
    
    % Parameters
    dim = dimensional_parameters();
    dim.Kd = 1; %1; %0; %1; %1000;
    dim.Ld = 100;
    scale = scalings(dim);
    pp = dimensionless_parameters(dim,scale);
    pp = isothermal_parameters(pp);
    % pp.B = 1000;
    pp.Kd = 10;
    
    % Initial condition
    phi00 = 0.2;
    Z0 = 0.5;
    phi0 = @(Z) phi00*exp(-Z/Z0); % max(0,0.2*(exp(-Z/0.5)-exp(-Z(end)/0.5))).*(Z<0.5); %"default";
    ss0 = isothermal_initial_cond(phi0,dd,pp);
    
    % Time span
    tspan = [0,1]; %[0,20];

elseif run_num == 2

    %% KAN_M PROMICE data
    % Discretisation
    dd = discretise_data_default();
    dd.plot = 0; %0;
    dd.gap = 10^3; %1000; %1;
    dd.dt = 1e-8; %3e-7; %dd.dt*1e-3; %dd.dt/200;

    % Parameters
    [dim,scale,pp] = real_forcing_parameters(dd);
    pp.Kd = 0; %1; %0.01;
    dim.Kd = pp.Kd*dim.Ld;
    
    % Initial conditions
    use_mean = 1;       % =1 to use mean forcing as intitial condition, =0 to use initial forcing
    ss0 = real_forcing_initial_cond(dd,pp,use_mean);
    ss0.S = 0.5*dd.dZ.^0; %0.5*dd.Z.^0;
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp); 
    
    % Time span
    tspan = [0,0.001]; %[0,1]; %[0,10];

elseif run_num == 3
    dd = discretise_isothermal();
    dd.plot = 0; %0;
    dd.gap = 1e2;
    dd.dt = 1e-7; %1e-8; %7e-5; %1e-6; %10*1e-5;

    dim = dimensional_parameters();
    Qsi = 200;
    Q0 = 50;
    dim.Qsi = @(t) Qsi;
    dim.Q0 = @(t) Q0;
    scale = scalings(dim);
    pp = dimensionless_parameters(dim,scale);
    pp.Kd = 0.1; %10; %0.1; %10;

    ss0 = steadily_melting_initial_cond(Qsi/scale.F,Q0/scale.F,pp,dd);  % Dimensionless Qsi and Q0 in here
    ss0.S = 1*dd.dZ.^0; %0.5*dd.Z.^0;
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp); 

    tspan = [0,1e-2]; %[0,0.1];
end

%% Timestepping
tic
ss = timestepping_unsat(ss0,pp,dd,tspan);
toc

%% Save data
save_time = inf;      % save the last save_time amount of time of data
save_ind = ss(end).t-[ss.t]<save_time;
ss_save = ss(save_ind);
tsave = [ss_save.t];

ff = struct();          % structure of forcing values (not relying on functions)
% ff.Qsi = scale.F*pp.Qsi(tsave);
% ff.Q0 = scale.F*pp.Q0(tsave);
% ff.nu = (scale.F/scale.T)*pp.nuhat(tsave);
% ff.a = a_func(tsave/scale.t);
% ff.alpha = alpha_func(tsave/scale.t);
% ff.r = r_func(tsave/scale.t);
save(["wc_evol_"+strrep(strrep(string(datetime),' ','_'),':','_')+".mat"],"ss_save","pp","dd","scale","dim","ff")

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
caxis([0,1])
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
