set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

addpath(genpath('./../..'))       % Add all subfolders to path
%%
% ss = ssh;

t = [ss.t];

Z0 = 0.5;
q0 = 0.1;
phi0 = 0.2;
% pp.Kd = 1;

tbegin = t(end)-60;

t = t(t>=tbegin);

cmap2 = colblind5();
Scol = cmap2(1,:);
phicol = 0*[1 1 1];
qcol = cmap2(4,:);
qdcol = cmap2(3,:);
Tcol = cmap2(5,:);
mintcol = cmap2(2,:);

%%
sequence_plot = figure(12);
set(gcf,'units','centimeters','paperpositionmode','auto','position',[11 2 7 7]); %,'renderer','painters');
clf;

for i=1:length(t)
    box on
    hold off
    % grid on
    p1 = plot(ss(i).phi*5,dd.Z,'color',phicol,'displayname','$\phi \times 5$'); %,'marker','x');
    hold on
    set(gca,'YDir','reverse')
    ylim([0,1])
    xlim([-1,1])
    p2 = plot(ss(i).S,dd.Z,'color',Scol,'displayname','$S$'); %,'marker','x');
    p5 = plot(ss(i).T,dd.Z,'color',Tcol,'displayname','$T$');
    p6 = plot(ss(i).mint/3,dd.Z,'color',mintcol,'DisplayName','$m_{\mathrm{int}}/3$');
    p3 = plot(ss(i).q*30,dd.Zedge,'color',qcol,'displayname','$q \times 30$');
    p4 = plot(ss(i).qd*20,dd.Z,'color',qdcol,'displayname','$q_\mathrm{d} \times 20$');

    title(['$t=',num2str(t(i) - tbegin,3),'$'])
    ylabel('$Z$','interpreter','latex');
    
    leg = legend('show','interpreter','latex','location','southwest'); %'position',[0.99,0.22,0.1,0.3]); %
    leg.ItemTokenSize = [10,1];
    
    xlabel(' ','HorizontalAlignment','left')

    if i == 1
        exportgraphics(gca,'snapshots.gif',Append=false)
    else
        exportgraphics(gca,'snapshots.gif',Append=true)
    end

end