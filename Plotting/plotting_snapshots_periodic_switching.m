set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

addpath(genpath('.'))       % Add all subfolders to path
%%
% ss = ssh;

t = [ss.t];

Z0 = 0.5;
q0 = 0.1;
phi0 = 0.2;
% pp.Kd = 1;

plot_mean_steady = 0;       % = 1 to plot the steadily melting solution of the mean forcing valuepp.Kd
plot_ana = 0; %1;       % 1 to plot the analytical solution too

tshift = 2;
period = 60;
tbegin = t(end) - tshift - period;
tstart = t(end) - tshift - period + 44.001; %26;
tend = t(end) - tshift - period + 62; %31; 
spacing = 1; %0.0001; %0.5; %2; 
% plot_times = t(end) - tshift - period + 0.001+ [22,22.4,22.8,23,24,25];%[tstart:spacing:tend]; % [40:1:45]; %[10,11,14,20,32]; %[[tstart:spacing:tend],t(end)];
% plot_times = [tstart:spacing:tend];
plot_times = t(end) - tshift - period + 0.001+ [49.9:0.1:50.4]; %[22,24,26,28,30]; %[44,45,45.5,46,48,51];
np = length(plot_times);
to_plot = zeros(np,1);
for i=1:np
    to_plot(i) = find(t>=plot_times(i),1);
end
%to_plot = length(ss)-[5,4,3,2,1,0];
plot_every = 1;     % plot every 'plot_every' of the plot_times times in the main figure - plot all in inset up to time 'inset_end'
inset_end = plot_times(end);     % final time to plot in the inset

cmap = colormap(parula(length(to_plot)+1));

cmap2 = colblind5();
Scol = cmap2(1,:);
phicol = 0*[1 1 1];
qcol = cmap2(4,:);
qdcol = cmap2(3,:);
Tcol = cmap2(5,:);
mintcol = cmap2(2,:);

%%
sequence_plot = figure(13);
set(gcf,'units','centimeters','paperpositionmode','auto','position',[11 2 18 5.5]); %,'renderer','painters');
clf;
ncols = 6; %ceil((length(to_plot)-1)/2);
tl = tiledlayout(1,ncols,'padding','compact','TileSpacing','compact');
% title(tl,['$\mathcal{K}_\mathrm{d}=',num2str(pp.Kd),'$'],'interpreter','latex')

mint_int = zeros(size(to_plot));
phi_diff_dt_int = zeros(size(to_plot));

for j=1:length(to_plot)
    mint_int(j) = trapz(dd.Z,ss(i).mint);
    phi_diff_dt_int(j) = trapz(dd.Z,(ss(i).phi-ss(i-1).phi)/(ss(i).t-ss(i-1).t));
ax(j) = nexttile;
box on
hold on
% grid on
set(gca,'YDir','reverse')
ylim([0,0.15])
xlim([-0.4,1]) %xlim([-0.1,1.2]) %([-0.1,2.3]) %xlim([0,1.1])
i = to_plot(j);
% p1 = plot(ss(i).phi*5,dd.Z,'color',phicol,'displayname','$\phi \times 5$','marker','x');
p1 = plot(ss(i).phi,dd.Z,'color',phicol,'displayname','$\phi$','marker','x');
p2 = plot(ss(i).S,dd.Z,'color',Scol,'displayname','$S$','marker','x');
p5 = plot(ss(i).T,dd.Z,'color',Tcol,'displayname','$T$');
p8 = plot(ss(i).H,dd.Z,'DisplayName','$H$');
p9 = plot(ss(i).W,dd.Z,'DisplayName','$W$');
% p6 = plot(ss(i).mint,dd.Z,'color',mintcol,'DisplayName','$m_{\mathrm{int}}$');
p6 = plot(ss(i).mint/10,dd.Z,'color',mintcol,'DisplayName','$m_{\mathrm{int}}/10$');
% p7 = plot((ss(i).phi-ss(i-1).phi)/(ss(i).t-ss(i-1).t),dd.Z,'Color',mintcol,'DisplayName','$(\phi-\phi_{\mathrm{old}})/ \mathrm{d} t$','LineStyle','--');
% p3 = plot(ss(i).q*30,dd.Zedge,'color',qcol,'displayname','$q \times 30$');
% p4 = plot(ss(i).qd*20,dd.Z,'color',qdcol,'displayname','$q_\mathrm{d} \times 20$');
if plot_ana==1
    [Sana,qana,qdana,phiana,Zb_func,Zt_func,Zw,ts] = isothermal_drainage_analytic(plot_times(j),dd.Z,pp.K,pp.Kd,phi0,Z0,q0);
    plot(Sana,dd.Z,'color',Scol,'linestyle','--','HandleVisibility','off');
    plot(10*qana,dd.Z,'color',qcol,'linestyle','--','HandleVisibility','off');
    plot(qdana,dd.Z,'color',qdcol,'linestyle','--','HandleVisibility','off');
end
title(['$t=',num2str(plot_times(j) - tbegin,3),'$'])
if j==1
ylabel('$Z$','interpreter','latex');
else
    yticklabels({})
end
if j==6
    leg = legend('show','interpreter','latex','location','northeast');
    leg.ItemTokenSize = [10,1];
end
if j==3
    xlabel(' ','HorizontalAlignment','left')
    % text(0.1,1.2,'Dimensionless variables','fontsize',11)
% elseif j==4
%     xlabel('variables','HorizontalAlignment','right')
end
end
linkaxes(ax)

% figure
% hold on
% plot(mint_int,'DisplayName','$m_{\mathrm{int}}$')
% plot(phi_diff_dt_int,'DisplayName','$(\phi-\phi_{\mathrm{old}})/ \mathrm{d} t$')
% legend('interpreter','latex')



%%
snap_plot = figure(14);
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
% print(sequence_plot,"periodic_switch_snapshots_"+strrep(strrep(string(datetime),' ','_'),':','_'),'-depsc','-r300')
% print(sequence_plot,"periodic_switch_snapshots_"+strrep(strrep(string(datetime),' ','_'),':','_'),'-dpng','-r300')
% print(sequence_plot,"periodic_switch_snapshots_"+strrep(strrep(string(datetime),' ','_'),':','_'),'-dsvg','-r300')