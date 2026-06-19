%%

set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

addpath(genpath('/scratch/woods/Dropbox/Oxford Mathematics DPhil/MATLAB/Weathering crust/Returning to WC Oct 2022'))
addpath(genpath('C:\Users\twood\Dropbox\Oxford Mathematics DPhil\MATLAB\Weathering crust\Returning to WC Oct 2022'))

addpath(genpath('/scratch/woods/Dropbox/Oxford Mathematics DPhil/MATLAB/colourmaps'))
addpath(genpath('C:\Users\twood\Dropbox\Oxford Mathematics DPhil\MATLAB\colourmaps'))

% load("wc_evol_24-Jun-2024_Kd0.mat")
load("wc_evol_02-Oct-2025_Kd0_B1e4.mat")
ssKd0 = ss_save;
ppKd0 = pp;
% load("wc_evol_24-Jun-2024_Kd2.mat")
load("wc_evol_02-Oct-2025_Kd2_B1e4.mat")
ssKd2 = ss_save;
ppKd2 = pp;
load("wc_evol_24-Jun-2024_Kd3.mat")
ssKd3 = ss_save;
ppKd3 = pp;
% load("wc_evol_24-Jun-2024_Kd10.mat")
load("wc_evol_02-Oct-2025_Kd10_B1e4.mat")
ssKd10 = ss_save;
ppKd10 = pp;

%%

ss = ssKd0;
pp = ppKd0;
Kd = pp.Kd;
t = [ss.t];

q0 = 0.1;
phi0 = 0.2;
Z0 = 0.5;

plot_mean_steady = 0;       % = 1 to plot the steadily melting solution of the mean forcing value
plot_ana = 1;       % 1 to plot the analytical solution too

tstart = 0.1; %t(1);
tend = t(end); %0.001; %t(end);
spacing = 0.1; %0.0001; %0.5; %2; 
plot_times = [tstart:spacing:tend]; %[[tstart:spacing:tend],t(end)];
np = length(plot_times);
to_plot = zeros(np,1);
for i=1:np
    to_plot(i) = find(t>=plot_times(i),1);
end
plot_every = 1;     % plot every 'plot_every' of the plot_times times in the main figure - plot all in inset up to time 'inset_end'
inset_end = plot_times(end);     % final time to plot in the inset

cmap = colblind8(); %colormap(parula(length(to_plot)+1));

cmap2 = colblind5();
Scol = cmap2(1,:);
phicol = 0.6*[1 1 1];
qcol = cmap2(4,:);
qdcol = cmap2(3,:);

Kd0_snap_plot = figure(11);
set(gcf,'units','centimeters','paperpositionmode','auto','position',[11 2 16 4.2],'renderer','painters');
clf;
ncols = 6; %ceil((length(to_plot)-1)/2);
tiledlayout(1,ncols,'TileSpacing','compact')

for j=1:length(to_plot)
nexttile
box on
hold on
% grid on
set(gca,'YDir','reverse')
ylim([0,1])
xlim([0,1]) %xlim([-0.1,1.2]) %([-0.1,2.3]) %xlim([0,1.1])
i = to_plot(j);
p1 = plot(ss(i).phi,dd.Z,'color',phicol,'displayname','$\phi$');
p2 = plot(ss(i).S,dd.Z,'color',Scol,'displayname','$S$');
p3 = plot(10*ss(i).q,dd.Zedge,'color',qcol,'displayname','$10q$');
if plot_ana==1
    [Sana,qana,qdana,phiana,Zb_func,Zt_func,Zw,ts] = isothermal_drainage_analytic(plot_times(j),dd.Z,pp.K,pp.Kd,phi0,Z0,q0);
    plot(Sana,dd.Z,'color',Scol,'linestyle','--','HandleVisibility','off');
    plot(10*qana,dd.Z,'color',qcol,'linestyle','--','HandleVisibility','off');
    % plot(qdana,dd.Z,'color',qdcol,'linestyle','--','HandleVisibility','off');
end
% p7 = plot([0,1],[ss(i).Zm,ss(i).Zm],'k--','displayname','$Z_\mathrm{m}$','handlevisibility','off');
% p11 = plot(ss(i).qd,dd.Z,'color',qdcol,'displayname','$10q_\mathrm{d}$');
title(['$t=',num2str(plot_times(j)),'$'])
% if or(j==1,j==6)
% ylabel('$Z$ (m)','interpreter','latex');
% end
if j==1
    leg = legend('show','interpreter','latex','location','southeast');
    leg.ItemTokenSize = [10,1];
    ylabel('$Z$')
    % text(-0.45,-0.15,'(b)')
    text(-0.8,0.5,['$\mathcal{K}_\mathrm{d}=',num2str(Kd),'$'],'rotation',90,'HorizontalAlignment','center')
end
if j>1
    yticks([])
end
end

%%
q0 = 0.1;

plot_mean_steady = 0;       % = 1 to plot the steadily melting solution of the mean forcing value
plot_ana = 1;       % 1 to plot the analytical solution too

tstart = t(1);
tend = 2; %0.001; %t(end);
spacing = 0.4; %0.0001; %0.5; %2; 
plot_times = [tstart:spacing:tend]; %[[tstart:spacing:tend],t(end)];
np = length(plot_times);

cmap = colblind8(); %colormap(parula(length(to_plot)+1));

cmap2 = colblind5();
Scol = cmap2(1,:);
phicol = 0.6*[1 1 1];
qcol = cmap2(4,:);
qdcol = cmap2(3,:);

Kd2_10_snap_plot = figure(12);
set(gcf,'units','centimeters','paperpositionmode','auto','position',[11 2 16 8],'renderer','painters');
clf;
ncols = 6; %ceil((length(to_plot)-1)/2);
tiledlayout(2,ncols,'padding','loose','TileSpacing','compact')

Kd_vals2 = [2,10];
for ii=1:length(Kd_vals2)
    Kd = Kd_vals2(ii);
    if Kd==2
        ss = ssKd2;
        pp = ppKd2;
    elseif Kd==10
        ss = ssKd10;
        pp = ppKd10;
    end
    t = [ss.t];
    to_plot = zeros(np,1);
    for i=1:np
        to_plot(i) = find(t>=plot_times(i),1);
    end

for j=1:length(to_plot)
nexttile
box on
hold on
% grid on
set(gca,'YDir','reverse')
ylim([0,1])
xlim([0,1]) %xlim([-0.1,1.2]) %([-0.1,2.3]) %xlim([0,1.1])
i = to_plot(j);
p1 = plot(ss(i).phi,dd.Z,'color',phicol,'displayname','$\phi$');
p2 = plot(ss(i).S,dd.Z,'color',Scol,'displayname','$S$');
p3 = plot(ss(i).q,dd.Zedge,'color',qcol,'displayname','$q$');
if plot_ana==1
    [Sana,qana,qdana,phiana,Zb_func,Zt_func,Zw,ts] = isothermal_drainage_analytic(plot_times(j),dd.Z,pp.K,pp.Kd,phi0,Z0,q0);
    plot(Sana,dd.Z,'color',Scol,'linestyle','--','HandleVisibility','off');
    plot(qana,dd.Z,'color',qcol,'linestyle','--','HandleVisibility','off');
    plot(qdana,dd.Z,'color',qdcol,'linestyle','--','HandleVisibility','off');
end
% p7 = plot([0,1],[ss(i).Zm,ss(i).Zm],'k--','displayname','$Z_\mathrm{m}$','handlevisibility','off');
p11 = plot(ss(i).qd,dd.Z,'color',qdcol,'displayname','$q_\mathrm{d}$');
if ii==1
title(['$t=',num2str(plot_times(j)),'$'])
xticklabels({})
end
% if or(j==1,j==6)
% ylabel('$Z$ (m)','interpreter','latex');
% end
if j==1
    if ii==length(Kd_vals2)
    leg = legend('show','interpreter','latex','location','southeast');
    leg.ItemTokenSize = [10,1];
    end
    ylabel('$Z$')
    % text(-0.45,-0.15,'(b)')
    text(-0.8,0.5,['$\mathcal{K}_\mathrm{d}=',num2str(Kd),'$'],'rotation',90,'HorizontalAlignment','center')
end
if j>1
    yticks([])
end
if and(ii == 1,j==length(to_plot))
    ss_inset = ss;
    qana_inset = qana;
end
end
end

% ss = ssKd2;
% pp = ppKd2;
% axes('position',[0.815,0.57,0.08,0.17])
axes('position',[0.815,0.57,0.08,0.17])
box on
hold on
p3 = plot(ss_inset(i).q,dd.Zedge,'color',qcol,'displayname','$q$');
plot(qana_inset,dd.Z,'color',qcol,'linestyle','--','HandleVisibility','off');
% p4 = plot(ss_inset(i).S,dd.Z,'color',Scol,'displayname','$S$');
set(gca,'ydir','reverse')
ylim([0,1])
xticks([])
yticks([])