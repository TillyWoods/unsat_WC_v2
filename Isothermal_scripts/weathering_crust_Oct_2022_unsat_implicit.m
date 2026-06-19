% Weathering crust code for when there is unsaturated flow and lateral
% drainage
% Feb 2024 (modified from October 2022)

% Used for both isothermal water transport test case and with real KAN_M
% forcing data

%%
addpath(genpath('/scratch/woods/Dropbox/Oxford Mathematics DPhil/MATLAB/Weathering crust/Returning to WC Oct 2022'))
addpath(genpath('C:\Users\twood\Dropbox\Oxford Mathematics DPhil\MATLAB\Weathering crust/Returning to WC Oct 2022'))
addpath(genpath('C:\Users\woods\Dropbox\Oxford Mathematics DPhil\MATLAB\Weathering crust/Returning to WC Oct 2022'))

%% Enter the run numner
% 1 for isothermal water transport
% 2 for KAN_M PROMICE data
% 3 for steadily melting solution with water flow
run_num = 1;

if run_num == 1
    %% Isothermal water transport
    % Discretisationpp.S
    dd = discretise_isothermal();

%     % Discretisation
% Zend = 10;         % size of domain
% Zfine = 3;          % size of fine domain near the surface
% dZedgeF = 0.01;     % grid size is fine part of domain
% dZedgeC = 1e-1;     % grid size in coarse part of domain
% NT = 10;            % number of grid cells in which to linearly transition between fine and coarse mesh
% Zedge = nonuniform_Zedge(Zend,Zfine,dZedgeF,dZedgeC,NT);
% 
% dd = discretise(Zedge);


    dd.plot = 0; %1;
    dd.dt = 0.5*0.5*1e-2;
    dd.gap = 10; %1e1; %1e2;
    
    % Parameters
    dim = dimensional_parameters();
    pp.Kd = 10; %2; %10;
    dim.Kd = pp.Kd*dim.Ld; %0; %1; %1000;
    % dim.Ld = 100;
    scale = scalings(dim);
    pp = dimensionless_parameters(dim,scale);
    pp = isothermal_parameters(pp);
    pp.tolS1 = 1e-10;

    pp.B = 1e4; %1e10;
    
    % Initial condition
    phi00 = 0.2;
    Z0 = 0.5;
    phi0 = @(Z) phi00*exp(-Z/Z0); % max(0,0.2*(exp(-Z/0.5)-exp(-Z(end)/0.5))).*(Z<0.5); %"default";
    % ss0 = ss04;
    ss0 = isothermal_initial_cond(phi0,dd,pp);
    
    % Time span
    tspan = [0,5]; %[0,0.41]; %[0,20];

elseif run_num == 2
    % Kd_vals = [0.003,0.1,0.01]; %[0.1,1,10,0.01,0.003];
    % for i=1:length(Kd_vals)
    %% KAN_M PROMICE data
    % Data
    load('wc_evol_KANM_SEB_2_vary_a.mat')
    ss_start = ss_save;
    t_sat = [ss_start.t_datetime];
    start_date = datetime('2008-09-02 14:00:00','Format','yyyy-MM-dd HH:mm:ss.SSS'); %datetime('2016-06-01 00:00:00','Format','yyyy-MM-dd HH:mm:ss.SSS'); % datetime('2016-06-01 00:00:00','Format','yyyy-MM-dd HH:mm:ss.SSS');
    end_date = start_date+years(10); %datetime('2016-09-15 00:00:00','Format','yyyy-MM-dd HH:mm:ss.SSS'); %datetime('2016-09-15 00:00:00','Format','yyyy-MM-dd HH:mm:ss.SSS');
    istart = find(t_sat>=start_date,1,'first');
    ss0 = ss_start(istart);
    ss0 = rmfield(ss0,'t_datetime');
    tend = days(end_date-start_date);

     % Discretisation
    dd = discretise_data_default();
    dd.plot = 0;
    dd.gap = 24; %2e2; %1e2; %10^3; %1000; %1;
    % dd.dt = 1e-6; %1e-6; %0.5e-5; %1e-8; %3e-7; %dd.dt*1e-3; %dd.dt/200;
    dd.dt = 1/24; % (1 hour) %1e-4; %0.25*1e-3;
    dd.steady_stop = 0; %1;
    dd.steady_tol = 1e-4; %1e-4;

    % Parameters
    tscale = "day";
    [dim,scale,pp] = real_forcing_parameters_tscale(dd,tscale,start_date);
    pp.Kd = 0; %1; %Kd_vals(i); %2; %1; %0.01;
    dim.Kd = pp.Kd*dim.Ld;

    pp.tolphi0 = 1e-5;
    
    % Initial conditions
    % use_mean = 1;       % =1 to use mean forcing as intitial condition, =0 to use initial forcing
    % ss0 = real_forcing_initial_cond(dd,pp,use_mean);
    % ss0.S = 0.01*dd.dZ.^0; %0.5*dd.Z.^0;
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp); 

    % ss0 = ss0_04;
    
    % Time span
    tspan = [0,tend]; %[0,1]; %[0,10];

elseif run_num == 3
    % Kd_vals = 3; %[10,3,2]; %[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,20,30,40,50,60,70,80,90]; %[0.1,100,3,4,5,6,7,8,9];
% 
    % for i=1:length(Kd_vals)
    %% Steadily melting with water flow
    dd = discretise_isothermal();
    dd.plot = 0;
    dd.gap = 1e2; %10; %1e2; %10; %1e2; %100;
    dd.dt = 0.25*1e-3; %1e-2; %1e-3; %0.25*1e-3; %0.25*1e-3; %1e-4; %7e-5; %1e-6; %10*1e-5;
    dd.steady_stop = 0;
    dd.steady_tol = 1e-4; %1e-4;

    dim = dimensional_parameters();
    Qsi = 100; %200; %0; %100; %100; %200;
    Q0 = 50; %-20; %20; %50; %100; 
    gap1 = 10; % 0; %10;
    gap2 = 20;
    dim.Qsi = @(t) 100*(t<gap1*dim.td); %200; %200*(t>=gap1*dim.td); %100*(t<gap1*dim.td); % Qsi;
    dim.Q0 = @(t) 50*(t<gap1*dim.td)-50*and(t>=gap1*dim.td,t<=(gap1+gap2)*dim.td)+100*(t>(gap1+gap2)*dim.td); % -20; %50*(t<gap1*dim.td)-50*and(t>=gap1*dim.td,t<=(gap1+gap2)*dim.td)+100*(t>(gap1+gap2)*dim.td); % Q0;
    % dim.Qsi = @(t) 100;
    % dim.Q0 = @(t) -100;
    scale = scalings(dim);
    pp = dimensionless_parameters(dim,scale);
    pp.Kd = 3; %0.003; %0.1; %Kd_vals(i); %1; %10; %0.1; %10;

    % pp.B = 1e4;

    pp.tolphi0 = 1e-4;

    ss0 = steadily_melting_initial_cond(Qsi/scale.F,Q0/scale.F,pp,dd);  % Dimensionless Qsi and Q0 in here
    % ss0.phi = 0*dd.Z.^0;
    ss0.S = 1*dd.Z.^0; %0.5*dd.Z.^0;
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp); 

    % ss0 = ss00;

    % ss0 = sspt1(end); %ssKd10_steady(end);
    tspan = [0,50]; %[0,1]; %[0,100]; %[0,0.1];
elseif run_num == 5
    %% Periodic switching - 30/05/2025
    % A periodic version of the switching experiment in the unsaturated
    % chapter of my thesis.
    % The forcing goes periodically between 3 sections of constant values.
    dd = discretise_isothermal();
    dd.plot = 0;
    dd.gap = 1e2; %10;
    dd.dt = 1e-2; %0.25*1e-3; %1e-2; %1e-3;
    dd.steady_stop = 0;
    dd.steady_tol = 1e-4; %1e-4;

    dim = dimensional_parameters();
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
    pp.Kd = 0; %0.003; %0.1; %Kd_vals(i); %1; %10; %0.1; %10;

    % pp.B = 1e4;

    pp.tolphi0 = 1e-4;

    ss0 = steadily_melting_initial_cond(Qsi/scale.F,Q0/scale.F,pp,dd);  % Dimensionless Qsi and Q0 in here
    % ss0.phi = 0*dd.Z.^0;
    ss0.S = 1*dd.Z.^0; %0.5*dd.Z.^0;
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp); 

    % ss0 = ss00;

    % ss0 = sspt1(end); %ssKd10_steady(end);
    tspan = [0,4*period]; %[0,1]; %[0,100]; %[0,0.1];

elseif run_num == 4
    %% Sinusoid
    % Idealised sinusoidal forcing

    % Kd_vals = 5; %[10,100,5,1];
    % for i=1:length(Kd_vals)

    dd = discretise_data_default();
    dd.gap = 1e3; %1e4; %10; %100; % 10 for 1e-3;
    dd.plot = 0;

    dd.dt = 0.25*1e-3; %1e-6; %1e-3;

    dim = dimensional_parameters();
    Qsi_A = 100; %100;
    Qsi_B = 0; %100;
    Q0_A = 50; %100;
    Q0_B = -50; %-50;
    Qsi_mean = (Qsi_A+Qsi_B)/2;
    Q0_mean = (Q0_A+Q0_B)/2;
    period = dim.ty;
    dim.Qsi = @(t) paused_sinusoid(t-period/4,Qsi_A,Qsi_B,period,0,0,0);
    dim.Q0 = @(t) paused_sinusoid(t-period/4,Q0_A,Q0_B,period,0,0,0);

    % dim.alpha = @(phi,S,t) dim.alpha0*(1+10*phi)./5;
    % dim.alpha = @(phi,S,t) dim.alpha0*(phi.*S+(1-phi));

    tscale = dim.td; %period; %dim.td; %dim.ty;
    scale = scalings_tscale(dim,tscale);
    pp = dimensionless_parameters(dim,scale);

    pp.Kd = Kd_vals(i); %100; %5;
    dd.periodic_stop = period/scale.t;
    dd.periodic_tol = 1e-3;

    ss0 = steadily_melting_initial_cond(Qsi_mean/scale.F,Q0_mean/scale.F,pp,dd);
    ss0.h = 0;
    % ss0.S = 0.5*dd.Z.^0;
    [ss0.H,ss0.W] = convert_T_to_H(ss0.T,ss0.phi,ss0.S,pp);

    % ss0 = ss_save(end);

    tspan = [0,2*365];
end

%% Timestepping
tic
ss = timestepping_unsat_implicit(ss0,pp,dd,tspan);
toc

%% Save data
save_time = inf;      % save the last save_time amount of time of data
save_ind = ss(end).t-[ss.t]<save_time;
% save_ind = 1:10:length(ss);
ss_save = ss(save_ind);
tsave = [ss_save.t];

%%

if run_num==2
for j = 1:length(ss_save)
    ss_save(j).t_datetime = start_date + days(ss_save(j).t);
end
end

% % ff = struct();          % structure of forcing values (not relying on functions)
% % ff.Qsi = scale.F*pp.Qsi(tsave);
% % ff.Q0 = scale.F*pp.Q0(tsave);
% % ff.nu = (scale.F/scale.T)*pp.nuhat(tsave);
% % ff.a = a_func(tsave/scale.t);
% % ff.alpha = alpha_func(tsave/scale.t);
% % ff.r = r_func(tsave/scale.t);
save(["wc_evol_"+strrep(strrep(string(datetime),' ','_'),':','_')+".mat"],"ss_save","pp","dd","scale","dim")
% save(["wc_evol_Kd"+num2str(Kd_vals(i))+"_"+strrep(strrep(string(datetime),' ','_'),':','_')+".mat"],"ss_save","pp","dd","scale","dim")
%     end
% end
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

function f = paused_sinusoid(t,val0,val1,period,pause0,pause1,pause2)
    % Creates a sinusoid that takes val0 at t=0, holds it there for time
    % pause0, then increases or
    % decreases smoothly to val1 as a sin/cos of period period until
    % t=period/2, then holds val1 for time pause1, then decreases or
    % increases back to val0 as a sinusoid, then holds val0 for time
    % pause2. Everything apart from the initial pasuse pause0 is then
    % repeated periodically.
    g = @(t) (val0+val1)/2+(val0-val1)/2*cos(2*pi*t/period);
    full_period = period + pause1 + pause2;
    t_mod = mod(t-pause0,full_period);
    log0 = and(t>=0,t<=pause0);
    log1 = and(t_mod>=0,t_mod <= period/2);
        f1 = g(t_mod);
    log2 = and(t_mod > period/2, t_mod < period/2 + pause1);
        f2 = val1*t.^0;
    log3 = and(t_mod >= period/2 + pause1,t_mod < period + pause1);
    log4 = t_mod >= period + pause1;
        f3 = g(t_mod-pause1);
        f4 = val0*t.^0;
    f = f1;
    f(log2) = f2(log2);
    f(log3) = f3(log3);
    f(log4) = f4(log4);
    f(log0) = f4(log0);
end
