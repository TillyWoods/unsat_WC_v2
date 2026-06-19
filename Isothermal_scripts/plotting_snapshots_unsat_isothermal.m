set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

addpath(genpath('/scratch/woods/Dropbox/Oxford Mathematics DPhil/MATLAB/Weathering crust/Returning to WC Oct 2022'))
addpath(genpath('C:\Users\twood\Dropbox\Oxford Mathematics DPhil\MATLAB\Weathering crust\Returning to WC Oct 2022'))

addpath(genpath('/scratch/woods/Dropbox/Oxford Mathematics DPhil/MATLAB/colourmaps'))
addpath(genpath('C:\Users\twood\Dropbox\Oxford Mathematics DPhil\MATLAB\colourmaps'))
%%
% ss = ssh;

t = [ss.t];

Z0 = 0.5;
q0 = 0.1;
phi0 = 0.2;
% pp.Kd = 1;

plot_mean_steady = 0;       % = 1 to plot the steadily melting solution of the mean forcing valuepp.Kd
plot_ana = 0; %1;       % 1 to plot the analytical solution too

tstart = t(1);
tend = t(end); %0.001; %t(end);
spacing = 0.2; %0.0001; %0.5; %2; 
plot_times = [tstart:spacing:tend]; %[[tstart:spacing:tend],t(end)];
np = length(plot_times);
to_plot = zeros(np,1);
for i=1:np
    to_plot(i) = find(t>=plot_times(i),1);
end
plot_every = 1;     % plot every 'plot_every' of the plot_times times in the main figure - plot all in inset up to time 'inset_end'
inset_end = plot_times(end);     % final time to plot in the inset

cmap = colormap(parula(length(to_plot)+1));

cmap2 = colblind5();
Scol = cmap2(1,:);
phicol = 0*[1 1 1];
qcol = cmap2(4,:);
qdcol = cmap2(3,:);

%%
sequence_plot = figure(13);
set(gcf,'units','centimeters','paperpositionmode','auto','position',[11 2 18 12]); %,'renderer','painters');
clf;
ncols = 5; %ceil((length(to_plot)-1)/2);
tl = tiledlayout(2,ncols,'padding','compact','TileSpacing','compact');
title(tl,['$\mathcal{K}_\mathrm{d}=',num2str(pp.Kd),'$'],'interpreter','latex')

for j=1:length(to_plot)
nexttile
box on
hold on
grid on
set(gca,'YDir','reverse')
ylim([0,1])
xlim([-0.1,1]) %xlim([-0.1,1.2]) %([-0.1,2.3]) %xlim([0,1.1])
i = to_plot(j);
p1 = plot(ss(i).phi,dd.Z,'color',phicol,'displayname','$\phi$');
p2 = plot(ss(i).S,dd.Z,'color',Scol,'displayname','$S$');
p3 = plot(ss(i).q*10,dd.Zedge,'color',qcol,'displayname','$10q$');
p4 = plot(ss(i).qd,dd.Z,'color',qdcol,'displayname','$q_\mathrm{d}$');
if plot_ana==1
    [Sana,qana,qdana,phiana,Zb_func,Zt_func,Zw,ts] = isothermal_drainage_analytic(plot_times(j),dd.Z,pp.K,pp.Kd,phi0,Z0,q0);
    plot(Sana,dd.Z,'color',Scol,'linestyle','--','HandleVisibility','off');
    plot(10*qana,dd.Z,'color',qcol,'linestyle','--','HandleVisibility','off');
    plot(qdana,dd.Z,'color',qdcol,'linestyle','--','HandleVisibility','off');
end
title(['$t=',num2str(plot_times(j)),'$'])
if or(j==1,j==6)
ylabel('$Z$','interpreter','latex');
end
if j==1
    leg = legend('show','interpreter','latex','location','southwest');
    leg.ItemTokenSize = [10,1];
end
end

%%
snap_plot = figure(13);
set(gcf,'units','centimeters','paperpositionmode','auto','position',[11 2 5 7]);
clf;
width = 1.2;

box on
hold on
set(gca,'YDir','reverse')
ylim([0,2.5])
xlim([0,1.1])

cmap = parula(length(to_plot)+1);

q0 = 0.1;
Sana = ((q0/pp.K)^(1/3))*(1./ss(1).phi).*exp(-pp.Kd*dd.Z/(pp.Ld*3));

for j=1:length(to_plot)
i = to_plot(j);
% p1 = plot(ss(i).phi,dd.Z,'k-') %,'displayname','$\phi$');
p2 = plot(ss(i).S,dd.Z,'color',cmap(j,:),'displayname',['$t=',num2str(t(i)),'$']); %b-','displayname','$S$');
% p3 = plot(10*ss(i).q,dd.Zedge,'g-') %,'displayname','$10q$');
end
% plot(ss(1).phi,dd.Z,'k-')
plot(Sana,dd.Z,'k--')
ylabel('$Z$ (m)','interpreter','latex');
xlabel('$S$')
% legend('show','interpreter','latex','location','southeast')

%%
% print(sequence_plot,"Figures/Just_made/unsat_snapshots_"+strrep(strrep(string(datetime),' ','_'),':','_'),'-depsc','-r300')
% print(sequence_plot,"Figures/Just_made/unsat_snapshots_"+strrep(strrep(string(datetime),' ','_'),':','_'),'-dpng','-r300')
% print(sequence_plot,"Figures/Just_made/unsat_snapshots_"+strrep(strrep(string(datetime),' ','_'),':','_'),'-dsvg','-r300')