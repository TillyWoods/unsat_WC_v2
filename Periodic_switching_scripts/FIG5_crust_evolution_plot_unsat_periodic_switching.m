%% Figure 5: periodic switching crust evolution plots
% November 2023
% Creates a figure with 4 panels from data produced by the weathering crust time evolution code (e.g. weathering_crust_Oct_2022). The panels are:
% •	Qsi and Q0 over time
% •	Temperature over depth and time
% •	Porosity over depth and time
% •	Surface melt rate, surface lowering rate, and crust thickness over time

addpath(genpath('.'))       % Add all subfolders to path

% Directory to save plots into
output_dir = ".";

set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

%% Load saved data, if desired (or just run after RUN_ME_periodic_switching_run.m)
% filename = "wc_periodic_switch_Kd_0_dZ_0-00125_dt_0-0003_14-Oct-2026_18_37_09.mat";  % no drainage
% filename = "wc_periodic_switch_Kd_20_dZ_0-00125_dt_2-5e-05_23-Oct-2026_22_49_27.mat";  % drainage
% load(filename)
ss = ss_save;

%%
set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

% Shift the period by a certain time for plotting
shift = 0; %5;

% Plot nondim or dim variables?
nondim = 1;     % 1 = plot nondim, 0 = plot dim

if nondim == 0
    scaleT = scale.T;
    scaleF = scale.F;
else
    scaleT = 1;
    scaleF = 1;
end

h = [ss.h];
phi = [ss.phi];
S = [ss.S];
T = [ss.T];
t1 = [ss.t];
hdot = [ss.hdot];
roff = [ss.roff];
msurf = [ss.msurf];
Zm = [ss.Zm];
H = [ss.H];
Zf = [ss.Zf];
% Zt = Zm - Zf;
Qsi = [ss.Qsi];
Q0 = [ss.Q0];

Z1 = ss(1).Z;
Zedge = ss(1).Zedge;

ztopT = -0;
ztopphi = -0;
zlimitT = -3; %-1.8; %-1.5; %-10; %-20; %-20;
zlimitphi = -1.6; %-1.8; %-1.5; %-1.5; %-0.3; %-5; %-13;
Ztop = 0;
Zlimit = 2;
tstart = t1(end)-60-shift; %t1(end)-60-shift;
tend = t1(end) - shift;
tstarti = find(t1>=tstart,1,'first');

red_lines = mod([20,40,60]+shift,60);    % where to put red horizontal lines, i.e. where the switching happens

fixed_surface = 0;      % 1 to use Z frame, 0 to use z frame
shift_z = 1;            % 1 to shift z so that h=0 at the start of the plot
shift_t = 1;            % 1 to shift t to be 0 at the start of the plot

if fixed_surface==0
    if shift_z==1
        indz = Z1<=-min(zlimitT,zlimitphi);
    else
        indz = Z1<=h(tstarti)-min(zlimitT,zlimitphi);
    end
else
    indz = Z1<=Zlimit;
end
indt = and(t1>=tstart,t1<=tend);
Z = Z1(indz);
t = t1(indt);
S = S(indz,indt);
phi = phi(indz,indt);
T = T(indz,indt);
H = H(indz,indt);
h = h(indt);
hdot = hdot(indt);
roff = roff(indt);
msurf = msurf(indt);
Zm = Zm(indt);
Zf = Zf(indt);
Qsi = Qsi(indt);
Q0 = Q0(indt);

[hgrid,Zgrid] = meshgrid(h,Z);
[tgrid,Zgrid] = meshgrid(t,Z);
zgrid = hgrid - Zgrid;

[hedgeg,Zedgeg] = meshgrid(h,Zedge);
[tedgeg,Zedgeg] = meshgrid(t,Zedge);
zedge = hedgeg - Zedgeg;

istart = find(t>=tstart,1);
iend = find(t<=tend,1,'last');

phi_max = min(max(max(phi(:,and(t>=tstart,t<=tend)))),1);
T_min = -3.5*(1-nondim)-0.4*nondim; %pp.Tinf*scaleT;
Q_max = max([scaleF*Qsi,scaleF*Q0])*1.1;
Q_min = min(0,min([scaleF*Qsi,scaleF*Q0]))*1.1;

if shift_z
    zgrid = zgrid-h(istart);
    h = h - h(istart);
    Zm = Zm + h(istart);
end

if shift_t
    % t = t(istart:iend) - tstart;
    t = t - tstart;
    tgrid = tgrid - tstart;
    % Zm = Zm(istart:iend);
    % h = h(istart:iend);
    % hdot = hdot(istart:iend);
    % msurf = msurf(istart:iend);
    % Zt = Zt(istart:iend);
    tmin = 0;
    tmax = tend-tstart;
    shift_t = tstart;
else
    % t = t(istart:iend);
    tmin = tstart;
    tmax = tend;
    shift_t = 0;
end

width = 1;

depth_grid = fixed_surface*Zgrid + (1-fixed_surface)*zgrid;

if scale.t == dim.td
    time_label = 'd';
elseif scale.t == dim.ty
    time_label = 'y';
end

xlinecol = 'r';
xlinestyle = '-';


if pp.Kd==0
    fignum = 20;
else
    fignum = 23;
end

Tphi_plot = figure(fignum); %(5);

clf
set(gcf,'units','centimeters','paperpositionmode','auto','position',[2 2 9 10*5/4]) %,'Renderer','painters'); %[2 2 0.8*9 0.8*10]); %
tiledlayout(5,1,'TileSpacing','compact','Padding','compact')

ax(3) = nexttile;
cols = parula(6);
col1 = cols(1,:);
col2 = cols(end-1,:);
box on
hold on
plot(t,scaleF*Qsi,'color',col2,'displayname','$Q_\mathrm{si}$')
plot(t,scaleF*Q0,'color',col1,'displayname','$Q_0$')
set(gca,'XTickLabel',{[]})
if nondim == 0
    ylabel('(W m$^{-2}$)')
end
leg = legend('show','interpreter','latex','position',[0.78,0.36+0.47,0.12,0.1-0.02]); %[0.78-0.05,0.36+0.47,0.2,0.1-0.02]); %
leg.ItemTokenSize = [15,10];
ylim([Q_min,Q_max])

colormap(parula(256));
snapshot_col = 0.5*[1 1 1];

if pp.Kd == 0
    title('(a) No drainage')
else
    title('(b) Drainage')
end

%%
ax(1) = nexttile;
box on
hold on
surf(tgrid,depth_grid,T*scaleT,'LineStyle','none')
view(0,90);
% xlim([tmin,tmax])
set(gca,'XTickLabel',{[]}) %set(gca,'XTick',[])
c1 = colorbar;
if nondim == 0
    Tlabel = '$T$ ($^{\mathrm{o}}$C)';
    Zlabel = '$Z$ (m)';
    zlabel = '$z$ (m)';
else
    Tlabel = '$T$';
    Zlabel = '$Z$';
    zlabel = '$z$';
end
ylabel(c1,Tlabel,'interpreter','latex')
if fixed_surface == 1
    plot3(t(Zm>0),Zm(Zm>0),max(T,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    plot3(t(Zf>0),Zf(Zf>0),max(T,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel(Zlabel,'interpreter','latex')
    set(gca,'YDir','reverse')
    % caxis([scale.T*pp.Tinf,0])
else
    if nondim == 0
        contourTvals = -flip([0:0.5:3]);
    else
        contourTvals = -flip([0:0.05:0.3]);
    end
    contour(tgrid,zgrid,T*scaleT,contourTvals,'color',0.5*[1 1 1])
    % plot3(t(Zm>0),h(Zm>0)-Zm(Zm>0),max(T,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    % plot3(t(Zf>0),h(Zf>0)-Zf(Zf>0),max(T,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel(zlabel,'interpreter','latex')
    plot3(t,h,max(T,[],'all')*t.^0,'k-','linewidth',width)
    caxis([T_min,0])
    contour(tgrid,depth_grid,T*scaleT,[0,0],'color','k')
    % contour(tgrid,depth_grid,phi,pp.tolphi0*[1 1],'color','k','linestyle','--')

    % contour3(tgrid,depth_grid,T*scaleT,flip(-[0:0.5:3]),'color',0.5*[1 1 1])
    
end
xline(red_lines,'color',xlinecol,'linestyle',xlinestyle)

%%
ax(2) = nexttile;
box on
hold on
surf(tgrid,depth_grid,phi,'LineStyle','none')
view(0,90);
% xlim([tmin,tmax])
set(gca,'XTickLabel',{[]}) %set(gca,'XTick',[])
c2 = colorbar;
caxis([0,phi_max])
ylabel(c2,'$\phi$','interpreter','latex')
if fixed_surface == 1
    plot3(t(Zm>0),Zm(Zm>0),max(phi,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    plot3(t(Zf>0),Zf(Zf>0),max(phi,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel(Zlabel,'interpreter','latex')
    set(gca,'YDir','reverse')
else
    % plot3(t(Zm>0),h(Zm>0)-Zm(Zm>0),max(phi,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    % plot3(t(Zf>0),h(Zf>0)-Zf(Zf>0),max(phi,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel(zlabel,'interpreter','latex')
    plot3(t,h,max(phi,[],'all')*t.^0,'k-','linewidth',width)
    contour3(tgrid,zgrid,phi,[0:0.05:1],'color',0.5*[1 1 1])

    contour(tgrid,depth_grid,phi,pp.tolphi0*[1 1],'color','k') %,'linestyle','--')
    % contour(tgrid,depth_grid,T*scaleT,[0,0],'color','k')
end
xline(red_lines,'color',xlinecol,'linestyle',xlinestyle)

ax(4) = nexttile;
box on
hold on
surf(tgrid,depth_grid,S,'LineStyle','none')
view(0,90);
set(gca,'XTickLabel',{[]}) %set(gca,'XTick',[])
% xlim([tmin,tmax])
% set(gca,'XTickLabel',{[]}) %set(gca,'XTick',[])
c2 = colorbar;
caxis([0,1])
ylabel(c2,'$S$','interpreter','latex')
if fixed_surface == 1
    plot3(t(Zm>0),Zm(Zm>0),max(S,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    plot3(t(Zf>0),Zf(Zf>0),max(S,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel(Zlabel,'interpreter','latex')
    set(gca,'YDir','reverse')
else
    % plot3(t(Zm>0),h(Zm>0)-Zm(Zm>0),max(S,[],'all')*t(Zm>0).^0,'k--','linewidth',width)
    % plot3(t(Zf>0),h(Zf>0)-Zf(Zf>0),max(S,[],'all')*t(Zf>0).^0,'k--','linewidth',width)
    ylabel(zlabel,'interpreter','latex')
    plot3(t,h,max(S,[],'all')*t.^0,'k-','linewidth',width)
    
    % contour3(tgrid,zgrid,S,[0:0.1:1],'color',0.5*[1 1 1])

    contour(tgrid,depth_grid,S,1e-10*[1 1],'color','k')
end
% xlabel(['$t$ (',time_label,')'])
xline(red_lines,'color',xlinecol,'linestyle',xlinestyle)

ax(5) = nexttile;
colKd = colblind8();
hdotcol = colKd(1,:);
msurfcol = colKd(2,:);
hold on
box on
plot(t,-hdot,'displayname','$-\dot{h}$','color',hdotcol)
% plot(t,roff,'r--','displayname','$r_\mathrm{off}$')
plot(t,msurf,'displayname','$m_\mathrm{surf}$','linestyle',':','color',msurfcol)
ylim([0,0.04])
if nondim == 0
    ylabel('(m d$^{-1}$)')
    xlabel(['$t$ (',time_label,')'])
else
    xlabel('$t$')
end
leg1 = legend('show','interpreter','latex','position',[0.78,0.1,0.23-0.08,0.07]);
leg1.ItemTokenSize = [10,1];

linkaxes(ax,'x')

if fixed_surface == 1
    ylim(ax(1),[Ztop,Zlimit])
    ylim(ax(2),[Ztop,Zlimit])
    ylim(ax(4),[Ztop,Zlimit])
else
    ylim(ax(1),[zlimitT,ztopT])
    ylim(ax(2),[zlimitphi,ztopphi])
    ylim(ax(4),[zlimitphi,ztopphi])
end

xlim([tmin,tmax])

%%
print(Tphi_plot,output_dir+"/evol_plot_Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-'),'-dpng','-r600')
print(Tphi_plot,output_dir+"/evol_plot_Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-'),'-depsc','-r600')
print(Tphi_plot,output_dir+"/evol_plot_Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-'),'-dsvg','-r600')