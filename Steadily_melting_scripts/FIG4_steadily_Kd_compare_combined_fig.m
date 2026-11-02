%% Steadily melting solutions plotting
% 21/06/2024
% Makes combined figure of how S, q, qd vary with Kd in steadily melting
% state, and shows how the saturated layer thickness varies, and that the
% total lateral drainage is the sum of the internal and surface melt

% Load saved steadily melting solutions for different lateral drainage Kd
load('ssKd_ssKdall_21-06-2024.mat')
%%

addpath(genpath('./..'))

set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

% Colours
n_cmap = 20;
cmap = sky(n_cmap);
cmap5 = colblind5();
cmap8 = colblind8();
mint_col = cmap8(5,:);
msurf_col = cmap8(3,:);
qint_col = cmap8(1,:);
cmap2 = colblind4();
col = cmap2(1,:);
colind = floor(Kd_vals)+1;
colind(end-1:end) = [n_cmap-1,n_cmap];

Sheight = 3;
Kdheight = 2;

smf = figure(1);
clf
set(gcf,'PaperPositionMode','auto','units','centimeters','position',[2 2 14 13],'Renderer','painters')
t = tiledlayout('flow','padding','compact','tilespacing','loose');

%% S
nexttile([Sheight,1])
box on
hold on
plot(ssKd(1).phi,ssKd(1).Z,'r-','displayname','$\phi$')
for i=1:length(Kd_vals)-1
    plot(ssKd(i).S,ssKd(i).Z,'displayname',['$\mathcal{K}_\mathrm{d}=',num2str(ssKd(i).Kd),'$'],'color',cmap(colind(i),:),'handlevisibility','off')
    if ssKd(i).Kd<100
    plot(ssKd(i).Sana,ssKd(i).Z,'linestyle','--','displayname',['$\mathcal{K}_\mathrm{d}=',num2str(ssKd(i).Kd),'$'],'color',cmap(colind(i),:),'handlevisibility','off')
    end
end
set(gca,'YDir','reverse')
xlabel('$S$')
ylabel('$Z$')
ylim([0,1.5])
annotation('arrow',[0.25,0.2],[0.8,0.6])
text(0.65,1.05,'$\mathcal{K}_\mathrm{d}$','interpreter','latex')
small_font = 8;
text(0.6,0.4,'$1$','interpreter','latex','fontsize',small_font)
text(0.3,0.82,'$10$','interpreter','latex','fontsize',small_font)
text(0.3,0.96,'$100$','interpreter','latex','fontsize',small_font)
text(0.88,0.3,'$0.1$','interpreter','latex','fontsize',small_font)
leg = legend('show','interpreter','latex');
leg.ItemTokenSize = [10,10];
title(['(a)',newline,'Saturation'])

%% q
nexttile([Sheight,1])
box on
hold on
for i=1:length(Kd_vals)-2
    p1 = plot(ssKd(i).q,ssKd(i).Zedge,'displayname',['$\mathcal{K}_\mathrm{d}=',num2str(ssKd(i).Kd),'$'],'color',cmap(colind(i),:),'handlevisibility','off');
    if ssKd(i).Kd<100
    p2 = plot(ssKd(i).qana,ssKd(i).Z,'linestyle','--','displayname',['$\mathcal{K}_\mathrm{d}=',num2str(ssKd(i).Kd),'$'],'color',cmap(colind(i),:),'handlevisibility','off');
    end
end
set(gca,'YDir','reverse')
xlabel('$q$')
ylabel('$Z$')
ylim([0,1.5])
annotation('arrow',[0.53,0.415],[0.81,0.905])
text(0.03,0.37,'$\mathcal{K}_\mathrm{d}$','interpreter','latex')
xlim([0,0.05])
title(['(b)',newline,'Downward flux'])

%% qd
nexttile([Sheight,1])
box on
hold on
for i=1:length(Kd_vals)-2
    p1 = plot(ssKd(i).qd,ssKd(i).Z,'displayname',['$\mathcal{K}_\mathrm{d}=',num2str(ssKd(i).Kd),'$'],'color',cmap(colind(i),:),'handlevisibility','off');
    if ssKd(i).Kd<100
    p2 = plot(ssKd(i).qdana,ssKd(i).Z,'linestyle','--','displayname',['$\mathcal{K}_\mathrm{d}=',num2str(ssKd(i).Kd),'$'],'color',cmap(colind(i),:),'handlevisibility','off');
    end
end
set(gca,'YDir','reverse')
xlabel('$q_\mathrm{d}$')
ylabel('$Z$')
ylim([0,1.5])
annotation('arrow',[0.74,0.89],[0.9,0.88])
text(0.18,0.14,'$\mathcal{K}_\mathrm{d}$','interpreter','latex')
leg = legend([p1,p2],{'Numerical','Semi-analytical'},'interpreter','latex','location','southeast');
leg.ItemTokenSize = [10,10];
title(['(c)',newline,'Lateral drainage'])


%% Saturated thickness as function of Kd
nexttile([Kdheight,2])
box on
hold on
plot([ssKdall.Kd],([ssKdall.Zm]-[ssKdall.Zd])./[ssKdall.Zm],'color',col,'displayname','Numerical')
plot([ssKdall.Kd],[ssKdall.dana]./[ssKdall.Zm],'color',col,'linestyle','--','displayname','Semi-analytical')
xlim([0,10])
xlabel('$\mathcal{K}_\mathrm{d}$')
ylabel('$d/Z_\mathrm{m}$')
title('(d)')
leg = legend('show','location','southwest','interpreter','latex');
leg.ItemTokenSize = [10,10];

%% Mass conservation bar chart
nexttile([Kdheight,1])
box on
hold on
Kdval = 5;
j = find([ssKdall.Kd]==Kdval,1,'first');
rangej = j;
b = bar([ssKdall(rangej).Kd],[ssKdall(rangej).mintint;ssKdall(rangej).msurf],'stacked');
b(1).FaceColor = mint_col;
b(2).FaceColor = msurf_col;
plot([ssKdall.Kd],[ssKdall.qdint],'color',qint_col,'linewidth',2)
ylim([0,0.05])
xlim([Kdval-0.6,Kdval+0.6])
title('(e)')
xticks([])
fsize = 9;
text(Kdval,0.005,['Subsurface',newline,'melt'],'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fsize)
text(Kdval,0.022,['Surface',newline,'melt'],'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fsize)
text(Kdval-0.18,0.042,['Total lateral',newline,'drainage'],'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fsize)
annotation('arrow',[0.8,0.74],[0.31,0.26])

%% Inset showing katurated thickness as function of Kd over larger Kd range
ax2 = axes('position',[0.4,0.22,0.2,0.1]);
box on
hold on
plot([ssKdall.Kd],([ssKdall.Zm]-[ssKdall.Zd])./[ssKdall.Zm],'color',col)
xlabel('$\mathcal{K}_\mathrm{d}$')
ylabel('$d/Z_\mathrm{m}$')
xlim([0,200])

%% Save fig
print(smf,"sm_fig_comb",'-dpng','-r600')
print(smf,"sm_fig_comb",'-depsc','-r600')
print(smf,"sm_fig_comb",'-dsvg','-r600')