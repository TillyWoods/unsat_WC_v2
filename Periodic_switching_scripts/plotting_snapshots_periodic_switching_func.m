function plotting_snapshots_periodic_switching_func(ss,pp,plot_times_raw,options)
% Function for plotting profiles (over Z) of various functions T, phi, S,
% and optionally also the water fluxes q and qd and the subsurface melt
% rate mint, at given snapshots in time.

arguments
    ss struct               % solution structure
    pp struct                % dimensionless parameters
    plot_times_raw double       % times of snapshots to plot
    options.fluxes logical = 0  % 1 to plot the water fluxes as well as T, phi, S
    options.mint logical = 0    % 1 to plot subsurface melt mint as well as T, phi, S
    options.output_dir string = "none"  % name of output directory, if you want the plot to be saved (as a png)
end

set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

% plot_times_raw used in paper:
% [22,22.4,22.8,23,24,25]
% [44,45,45.5,46,48,51]

base_name = "Kd_"+num2str(pp.Kd,2)+"_t_"+strrep(num2str(plot_times_raw(1))+"_dZ_"+strrep(num2str(ss(1).Zedge(2)-ss(1).Zedge(1)),'.','-'),'.','-');

addpath(genpath('.'))       % Add all subfolders to path

t = [ss.t];
Z = ss(1).Z;
Zedge = ss(1).Zedge;

tshift = 0;
period = 60;
tbegin = t(end) - tshift - period;
plot_times = t(end) - tshift - period + plot_times_raw; %+0.001;

np = length(plot_times);
to_plot = zeros(np,1);
for i=1:np
    to_plot(i) = find(t>=plot_times(i),1);
end
plot_times = t(to_plot);

cmap2 = colblind5();
Scol = cmap2(1,:);
phicol = 0*[1 1 1];
qcol = cmap2(4,:);
qdcol = cmap2(3,:);
Tcol = cmap2(5,:);
mintcol = cmap2(2,:);

%%
sequence_plot = figure(12);
set(gcf,'units','centimeters','paperpositionmode','auto','position',[11 2 18 5.5]); %,'renderer','painters');
clf;
ncols = 6; %ceil((length(to_plot)-1)/2);
tl = tiledlayout(1,ncols,'padding','compact','TileSpacing','compact');
% title(tl,['$\mathcal{K}_\mathrm{d}=',num2str(pp.Kd),',  \Delta Z =',num2str(Zedge(2)-Zedge(1)),'$'],'interpreter','latex')

mint_int = zeros(size(to_plot));
phi_diff_dt_int = zeros(size(to_plot));

for j=1:length(to_plot)
    mint_int(j) = trapz(Z,ss(i).mint);
    phi_diff_dt_int(j) = trapz(Z,(ss(i).phi-ss(i-1).phi)/(ss(i).t-ss(i-1).t));
ax(j) = nexttile;
box on
hold on
% grid on
set(gca,'YDir','reverse')
ylim([0,1])
xlim([-0.4,1]) %xlim([-0.1,1.2]) %([-0.1,2.3]) %xlim([0,1.1])
i = to_plot(j);
p1 = plot(ss(i).phi*5,Z,'color',phicol,'displayname','$\phi \times 5$'); %,'marker','x');
% p1 = plot(ss(i).phi,Z,'color',phicol,'displayname','$\phi$','marker','x');
p2 = plot(ss(i).S,Z,'color',Scol,'displayname','$S$'); %,'marker','x');
p5 = plot(ss(i).T,Z,'color',Tcol,'displayname','$T$');
% p8 = plot(ss(i).H,Z,'DisplayName','$H$');
% p9 = plot(ss(i).W,Z,'DisplayName','$W$');

if options.mint == 1
    % p6 = plot(ss(i).mint,Z,'color',mintcol,'DisplayName','$m_{\mathrm{int}}$');
    p6 = plot(ss(i).mint/5,Z,'color',mintcol,'DisplayName','$m_{\mathrm{int}}/5$');
    % p7 = plot((ss(i).phi-ss(i-1).phi)/(ss(i).t-ss(i-1).t),Z,'Color',mintcol,'DisplayName','$(\phi-\phi_{\mathrm{old}})/ \mathrm{d} t$','LineStyle','--');
end
if options.fluxes == 1 
    p3 = plot(ss(i).q*30,Zedge,'color',qcol,'displayname','$q \times 30$');
    p4 = plot(ss(i).qd*20,Z,'color',qdcol,'displayname','$q_\mathrm{d} \times 20$');
end

title(['$t=',num2str(plot_times(j) - tbegin,3),'$'])
if j==1
ylabel('$Z$','interpreter','latex');
else
    yticklabels({})
end
if j==6
    leg = legend('show','interpreter','latex','position',[0.89,0.22,0.1,0.3]); %'location','northeast');
    leg.ItemTokenSize = [10,1];
end
if j==3
    xlabel(' ','HorizontalAlignment','left')
end
end
linkaxes(ax)

if options.output_dir ~= "none"
    print(sequence_plot,options.output_dir+"/periodic_switch_snapshots_"+base_name,'-depsc','-r300')
    print(sequence_plot,options.output_dir+"/periodic_switch_snapshots_"+base_name,'-dpng','-r600')
    % print(sequence_plot,options.output_dir+"/periodic_switch_snapshots_"+base_name,'-dsvg','-r300')
end

end