%% Subsurface vs. surface melt and runoff for periodic switching
% 30/05/2025
% Modified 15/10/2026 to make a single plot, i.e. not comparing drainage
% and no drainage

addpath(genpath('.'))       % Add all subfolders to path

set(0,'DefaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultColorbarTickLabelInterpreter','latex')
set(0,'DefaultLineLineWidth',1)

% Directory to save plots into
output_dir = "OutputPlots_16-10-2026_tols"; %"OutputPlots_14-10-2026"; %

shift = 0; %20;%5;
trange = @(t) and(t<=t(end)-shift,t>t(end)-60-shift); %and(t<=t(end)-60-shift,t>t(end)-120-shift); %

t = [ss.t];
ss = ss(trange(t));

% Calculate subsurface vs. surface melt and runoff
ssmr = surf_int_percentages(ss);

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

bars = figure(15);
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[2 2 7 6]) %6])
clf
% box on
hold on
i=1;
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
if pp.Kd > 0
text(2*i+pshift-barshift,(ssmr.qd_cum_end+ssmr.freeze_cum_end)-vs,[num2str(ssmr.qdp*100,'%.0f'),'\%'],'HorizontalAlignment','right','fontsize',pfs,'color','w')
end
text(2*i+pshift-barshift,(ssmr.qd_cum_end+ssmr.freeze_cum_end+ssmr.roff_cum_end)-vs,[num2str(ssmr.roffp*100,'%.0f'),'\%'],'HorizontalAlignment','right','fontsize',pfs,'color','k')

xticks([1.5]) %,3.5])
% xticklabels({'$\mathcal{K}_{\mathrm{d}}=0$','$\mathcal{K}_{\mathrm{d}}=3$'})
if pp.Kd==0
    xticklabels({'\textbf{No drainage}'})
else
    xticklabels({'\textbf{Drainage}'})
end
% yticks([])
ylabel('Meltwater equivalent') % (m)')
xlim([0.5,2.4]) %4.5])
ylim([-0.05,1.08])


text(1+barshift,1.06,['MELT'],'HorizontalAlignment','center','color','k','fontsize',fs-1)
text(1+barshift,0.6,['Surface'],'HorizontalAlignment','center','color','k','fontsize',fs)
text(1+barshift-2*barshift,0.04,['Subsurface'],'HorizontalAlignment','center','color','k','fontsize',fs)
%text(3+barshift,1.06,['MELT'],'HorizontalAlignment','center','color','k','fontsize',fs-1)
if pp.Kd>0
    text(2-barshift,0.28,['Subsurface',newline,'lateral',newline,'drainage'],'HorizontalAlignment','center','color','w','fontsize',fs)
end
text(2-barshift,0.75,['Surface',newline,'runoff'],'HorizontalAlignment','center','color','k','fontsize',fs)
text(2-barshift,1.06,['REMOVAL'],'HorizontalAlignment','center','color','k','fontsize',fs-1)
%text(4-barshift,1.06,['REMOVAL'],'HorizontalAlignment','center','color','k','fontsize',fs-1)
% text(2,78-3,['Surface',newline,'runoff'],'HorizontalAlignment','center','color','k','fontsize',fs)
% text(2,54-3,['Refreezing'],'HorizontalAlignment','center','color','k','fontsize',fs)
text(2-barshift-0*barshift,0.02,['Refreezing'],'HorizontalAlignment','center','color','w','fontsize',fs)
% text(2.4,-0.04,['Refreezing'],'HorizontalAlignment','center','color','k','fontsize',fs)

%% Checking about the numerical mass loss
lower = ss(1).h-ss(end).h;
roff_tot = ssmr.roff_cum_end+ssmr.qd_cum_end;

% print(bars,output_dir+"/bars_Kd_"+num2str(pp.Kd,2)+"_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-'),'-dpng','-r600')
% print(bars,output_dir+"/bars_Kd_"+num2str(pp.Kd,2)+"_NO_TRAPZ_dZ_"+strrep(num2str(dd.dZ(2)),'.','-')+"_dt_"+strrep(num2str(dd.dt),'.','-'),'-dpng','-r600')
%print(bars,output_dir+"/bars_Kd_"+num2str(pp.Kd,2)+"_tolphi0_"+strrep(num2str(pp.tolphi0),'.','-')+"_tolS1_"+strrep(num2str(pp.tolS1),'.','-'),'-dpng','-r600')


function ssmr = surf_int_percentages(ss)
    ssmr.Z = ss(1).Z;
    ssmr.t = [ss.t]-ss(1).t;
    Zedge = ss(1).Zedge;

    % if isfield(ss,"qd_cum")
    %     ssmr.qd_cum = trapz(ssmr.Z,[ss.qd_cum]-ss(1).qd_cum);
    % else
    %     ssmr.qd_cum = zeros(size(ssmr.t));
    % end
    % ssmr.mint_cum = trapz(ssmr.Z,[ss.mint_cum]-ss(1).mint_cum);
    % ssmr.msurf_cum = [ss.msurf_cum]-ss(1).msurf_cum;
    % ssmr.roff_cum = [ss.roff_cum]-ss(1).roff_cum;
    % ssmr.freeze_cum = -trapz(ssmr.Z,[ss.freeze_cum]-ss(1).freeze_cum);

    top_id = 2;
    
    if isfield(ss,"qd_cum")
        qd_cum = [ss.qd_cum];
        ssmr.qd_cum = sum((qd_cum(top_id:end,:)-qd_cum(top_id:end,1)).*diff(Zedge(top_id:end))); %trapz(ssmr.Z,[ss.qd_cum]-ss(1).qd_cum);
    else
        ssmr.qd_cum = zeros(size(ssmr.t));
    end

    meltin_cum = [ss.mint_cum];
    freeze_cum = [ss.freeze_cum];
    ssmr.mint_cum = sum((meltin_cum(top_id:end,:)-meltin_cum(top_id:end,1)).*diff(Zedge(top_id:end))); %trapz(ssmr.Z,[ss.mint_cum]-ss(1).mint_cum);
    ssmr.msurf_cum = [ss.msurf_cum]-ss(1).msurf_cum;
    ssmr.roff_cum = [ss.roff_cum]-ss(1).roff_cum;
    ssmr.freeze_cum = -sum((freeze_cum(top_id:end,:)-freeze_cum(top_id:end,1)).*diff(Zedge(top_id:end))); %-trapz(ssmr.Z,[ss.freeze_cum]-ss(1).freeze_cum);
    
    % ssmr.mint_cum = sum(([ss.mint_cum]-ss(1).mint_cum).*diff([ss(1).Zedge])); %trapz(ssmr.Z,[ss.mint_cum]-ss(1).mint_cum);
    % ssmr.msurf_cum = [ss.msurf_cum]-ss(1).msurf_cum;
    % ssmr.roff_cum = [ss.roff_cum]-ss(1).roff_cum;
    % ssmr.freeze_cum = -sum(([ss.freeze_cum]-ss(1).freeze_cum).*diff([ss(1).Zedge])); %-trapz(ssmr.Z,[ss.freeze_cum]-ss(1).freeze_cum);
    
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