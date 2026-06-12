%% Internal vs. surface melt and runoff for periodic switching
% 30/05/2025

addpath(genpath('.'))       % Add all subfolders to path

set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

shift = 0; %20;%5;
trange = @(t) and(t<=t(end)-shift,t>t(end)-60-shift);

% Drainage (Kd=20)
load("../wc_periodic_switch_Kd_20_11-Jun-2026_dZ_0-0025_dt_halved.mat")
ssd = ss_save;
t = [ssd.t];
ssd = ssd(trange(t));

% No drainage (Kd=0)
load("../wc_periodic_switch_Kd_0_11-Jun-2026_dZ_0-0025_dt_halved.mat")

ssnd = ss_save;
t = [ssnd.t];
% ssnd = ssnd(t>t(end)-60);
ssnd = ssnd(trange(t));

% Calculate internal vs. surface melt and runoff
ssmrd = surf_int_percentages(ssd);
ssmrnd = surf_int_percentages(ssnd);

%%
cmap8 = colblind8();
mint_col = cmap8(5,:);
msurf_col = cmap8(3,:);
qd_col = cmap8(1,:);
freeze_col = cmap8(6,:);
roff_col = cmap8(4,:);
loss_col = 0.8*[1 1 1];

%%
fs = 10;
pfs = 8; %9;

pshift = 0.4; %0.19;

vs = 0.04;

barshift = 0.05;

barwidth = 0.84;

figure(12)
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[2 2 13 6]) %6])
clf
% box on
hold on
% ssmr = ssmrd;
% b = bar(0,[ssmr.mintp,ssmr.msurfp]*100,'stacked');
% b(1).FaceColor = mint_col;
% b(2).FaceColor = msurf_col;
% text(0+pshift,ssmr.mintp*100-5,[num2str(ssmr.mintp*100,'%.0f'),'\%'],'fontsize',pfs,'color','k')
% text(0+pshift,(ssmr.mintp+ssmr.msurfp)*100-5,[num2str(ssmr.msurfp*100,'%.0f'),'\%'],'fontsize',pfs,'color','k')
for i = 1:2
    if i==1
        ssmr = ssmrnd;
    elseif i==2
        ssmr = ssmrd;
    end
    b = bar(2*i-1+barshift,[ssmr.mint_cum_end,ssmr.msurf_cum_end],barwidth,'stacked');
    b(1).FaceColor = mint_col;
    b(2).FaceColor = msurf_col;
    b = bar(2*i-barshift,[ssmr.freeze_cum_end,ssmr.qd_cum_end,ssmr.roff_cum_end,ssmr.mass_loss],barwidth,'stacked');
    b(1).FaceColor = freeze_col;
    b(2).FaceColor = qd_col;
    b(3).FaceColor = roff_col;
    b(4).FaceColor = loss_col;
    text(2*i+pshift+2.5*barshift,ssmr.freeze_cum_end-0.5*vs,[num2str(ssmr.freezep*100,'%.0f'),'\%'],'HorizontalAlignment','right','fontsize',pfs,'color','k')
    text(2*i-1+pshift+barshift,ssmr.mint_cum_end-vs,[num2str(ssmr.mintp*100,'%.0f'),'\%'],'HorizontalAlignment','right','fontsize',pfs,'color','k')
    text(2*i-1+pshift+barshift,(ssmr.mint_cum_end+ssmr.msurf_cum_end)-vs,[num2str(ssmr.msurfp*100,'%.0f'),'\%'],'HorizontalAlignment','right','fontsize',pfs,'color','k')
    if i>1
    text(2*i+pshift-barshift,(ssmr.qd_cum_end+ssmr.freeze_cum_end)-vs,[num2str(ssmr.qdp*100,'%.0f'),'\%'],'HorizontalAlignment','right','fontsize',pfs,'color','w')
    end
    text(2*i+pshift-barshift,(ssmr.qd_cum_end+ssmr.freeze_cum_end+ssmr.roff_cum_end)-vs,[num2str(ssmr.roffp*100,'%.0f'),'\%'],'HorizontalAlignment','right','fontsize',pfs,'color','k')
end
xticks([1.5,3.5])
% xticklabels({'$\mathcal{K}_{\mathrm{d}}=0$','$\mathcal{K}_{\mathrm{d}}=3$'})
xticklabels({'\textbf{No drainage}','\textbf{Drainage}'})
% yticks([])
ylabel('Meltwater equivalent') % (m)')
xlim([0.5,4.5])


text(1+barshift,1.06,['MELT'],'HorizontalAlignment','center','color','k','fontsize',fs-1)
text(1+barshift,0.6,['Surface'],'HorizontalAlignment','center','color','k','fontsize',fs)
text(1+barshift-2*barshift,0.04,['Internal'],'HorizontalAlignment','center','color','k','fontsize',fs)
text(3+barshift,1.06,['MELT'],'HorizontalAlignment','center','color','k','fontsize',fs-1)
text(4-barshift,0.28,['Subsurface',newline,'lateral',newline,'drainage'],'HorizontalAlignment','center','color','w','fontsize',fs)
text(2-barshift,1.06,['REMOVAL'],'HorizontalAlignment','center','color','k','fontsize',fs-1)
text(4-barshift,1.06,['REMOVAL'],'HorizontalAlignment','center','color','k','fontsize',fs-1)
% text(2,78-3,['Surface',newline,'runoff'],'HorizontalAlignment','center','color','k','fontsize',fs)
% text(2,54-3,['Refreezing'],'HorizontalAlignment','center','color','k','fontsize',fs)
text(2-barshift,0.55,['Surface',newline,'runoff'],'HorizontalAlignment','center','color','k','fontsize',fs)
text(2-barshift-0*barshift,0.02,['Refreezing'],'HorizontalAlignment','center','color','w','fontsize',fs)

% text(2.4,-0.04,['Refreezing'],'HorizontalAlignment','center','color','k','fontsize',fs)

%% Checking about the numerical mass loss
lowernd = ssnd(1).h-ssnd(end).h;
roff_totnd = ssmrnd.roff_cum_end+ssmrnd.qd_cum_end;



function ssmr = surf_int_percentages(ss)
    ssmr.Z = ss(1).Z;
    ssmr.t = [ss.t]-ss(1).t;
    if isfield(ss,"qd_cum")
        ssmr.qd_cum = trapz(ssmr.Z,[ss.qd_cum]-ss(1).qd_cum);
    else
        ssmr.qd_cum = zeros(size(ssmr.t));
    end
    ssmr.mint_cum = trapz(ssmr.Z,[ss.mint_cum]-ss(1).mint_cum);
    ssmr.msurf_cum = [ss.msurf_cum]-ss(1).msurf_cum;
    ssmr.roff_cum = [ss.roff_cum]-ss(1).roff_cum;
    ssmr.freeze_cum = -trapz(ssmr.Z,[ss.freeze_cum]-ss(1).freeze_cum);
    
    if isfield(ss,"qd_cum")
        ssmr.qd_cum_end = ssmr.qd_cum(end);
    else
        ssmr.qd_cum_end = 0;
    end
    ssmr.mint_cum_end = ssmr.mint_cum(end);
    ssmr.msurf_cum_end = ssmr.msurf_cum(end);
    ssmr.roff_cum_end = ssmr.roff_cum(end);
    ssmr.freeze_cum_end = ssmr.freeze_cum(end);

    ssmr.melt_tot = ssmr.mint_cum_end + ssmr.msurf_cum_end;
    ssmr.loss_tot = ssmr.roff_cum_end + ssmr.qd_cum_end + ssmr.freeze_cum_end;
    ssmr.mass_loss = ssmr.melt_tot - ssmr.loss_tot;
    ssmr.roff_tot = ssmr.roff_cum_end+ssmr.qd_cum_end;

    ssmr.msurfp = ssmr.msurf_cum_end/ssmr.melt_tot;
    ssmr.mintp = ssmr.mint_cum_end/ssmr.melt_tot;
    ssmr.roffp = ssmr.roff_cum_end/ssmr.melt_tot;
    ssmr.freezep = ssmr.freeze_cum_end/ssmr.melt_tot;
    if isfield(ss,"qd_cum")
        ssmr.qdp = ssmr.qd_cum_end/ssmr.melt_tot;
    else
        ssmr.qdp = 0;
    end
    ssmr.mass_lossp = ssmr.mass_loss/ssmr.melt_tot;
end