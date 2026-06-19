% 24/06/2024
% Making figure comparing evolution of isothermal water transport
% saturation with semi-analytical and numerical solutions.

% addpath(genpath('/scratch/woods/Dropbox/Oxford Mathematics DPhil/MATLAB/matlab-epsclean-master'))
% addpath(genpath('C:\Users\twood\Dropbox\Oxford Mathematics DPhil\MATLAB\matlab-epsclean-master'))

%%
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

%% Colour plots of solutions combined
t = linspace(0,5,500);
Z = linspace(0,2.5,250);

phi0 = 0.2;

Z0 = 0.5;
q0 = 0.1;
K = 100;

Kdmax = 3/Z0;

Kd_vals = [0,2,10];

dratio = 1; %1.5; %2;

Kdevol_fig = figure(13);
clf
set(gcf,'units','centimeters','paperpositionmode','auto','position',[10 10 18 8]) %,'render','painters')
% tiledlayout(1,3,'padding','compact')

h = 0.4;
tend_max = 4;

y0num = 0.1;
y0ana = y0num + 1.1*h;
gaps = 0.02;
widths = [0.1,0.26,0.35+0.09]*1.08;
widths = [0.1,0.26,0.26+0.073]*1.25;
widths = [0.1,0.26*tend_max/2.2,(0.26+0.07)*tend_max/2.2]*1.25*0.48*1.2;
x0 = [2*gaps,3*gaps+widths(1),4*gaps+widths(2)+widths(1)]+gaps+0.03;
w = [widths]; %0.2,0.3,0.4];

tends = zeros(length(Kd_vals));

% clear ax;

for i=1:length(Kd_vals)
    
    % Analytical
    Kd = Kd_vals(i);
[S,q,qd,phi,Zb_func,Zt_func,Zw,ts]=isothermal_drainage_analytic(t,Z,K,Kd,phi0,Z0,q0);

% w = h*Zend/ts/dratio;
% ax(i) =
axes('position',[x0(i),y0ana,w(i),h],'units','normalized');
% nexttile
% x0 = x0+w;
box on
hold on

sat_surf = find(S(1,:)==1,1,'first');
sat_surf(isempty(sat_surf)) = length(t);

Zend = 2; %Z(end);
tend = min(t(sat_surf),tend_max); %2.2);
tends(i) = tend;

[tgrid,Zgrid] = meshgrid(t,Z);
surf(tgrid,Zgrid,S,'LineStyle','none')
view(0,90);
line(t(t<ts),Zw(t(t<ts)),2*t(t<ts).^0,'linewidth',1,'color','red','linestyle','-')
if ~isinf(ts)
line(t(t>ts),Zt_func(t(t>ts)),2*t(t>ts).^0,'linewidth',1,'color','k','linestyle','-')
line(t(t>ts),Zb_func(t(t>ts)),2*t(t>ts).^0,'linewidth',1,'color',['k'],'linestyle','-')
end
daspect([1,dratio,1])
xlim([min(t),tend])
ylim([0,Zend])
% xlabel('$t$','interpreter','latex')
xticklabels({})
if i==1
ylabel('$Z$','interpreter','latex')
text(-0.8,Zend/2,'Semi-analytical','Rotation',90,'HorizontalAlignment','center')
else
    yticks([])
end
title(['$\mathcal{K}_\mathrm{d}=',num2str(Kd),'$'],'interpreter','latex')
set(gca,'YDir','reverse')
if i==length(Kd_vals)
c1 = colorbar;
caxis([0,1])
ylabel(c1,'$S$','interpreter','latex')
end

% Numerical
    if Kd==0
        ss = ssKd0;
    elseif Kd ==2
        ss = ssKd2;
    elseif Kd==3
        ss = ssKd3;
    elseif Kd==10
        ss = ssKd10;
    end
    t1 = [ss.t];
    Z1 = [ss(end).Z];
    S1 = [ss.S];

% w = h*Zend/ts/dratio;
% ax(i+3) =
axes('position',[x0(i),y0num,w(i),h],'units','normalized');
box on
hold on

tnum = t1(t1<=tend*1.1);
Znum = Z1(Z1<=Zend);
Snum = S1(Z1<=Zend,t1<=tend*1.1);

[tgridnum,Zgridnum] = meshgrid(tnum,Znum);
surf(tgridnum,Zgridnum,Snum,'LineStyle','none')
view(0,90);
line(t(t<ts),Zw(t(t<ts)),2*t(t<ts).^0,'linewidth',1,'color','red','linestyle','-')
if ~isinf(ts)
line(t(t>ts),Zt_func(t(t>ts)),2*t(t>ts).^0,'linewidth',1,'color','k','linestyle','-')
line(t(t>ts),Zb_func(t(t>ts)),2*t(t>ts).^0,'linewidth',1,'color',['k'],'linestyle','-')
end
daspect([1,dratio,1])
xlim([min(t),tend])
ylim([0,Zend])
xlabel('$t$','interpreter','latex')
if i==1
ylabel('$Z$','interpreter','latex')
text(-0.8,Zend/2,'Numerical','Rotation',90,'HorizontalAlignment','center')
else
    yticks([])
end
% title(['$\mathcal{K}_\mathrm{d}=',num2str(Kd),'$'],'interpreter','latex')
set(gca,'YDir','reverse')
if i==length(Kd_vals)
c1 = colorbar;
caxis([0,1])
ylabel(c1,'$S$','interpreter','latex')
end

end

%%
% exportgraphics(Kdevol_fig,'iso_ana_num_comp_evol2.eps','ContentType','vector','BackgroundColor','none','Resolution',600)
% %%
% exportgraphics(Kdevol_fig,'iso_ana_num_comp_evol2.png','Resolution',600)
% %%
% exportgraphics(Kdevol_fig,'iso_ana_num_comp_evol2.svg','ContentType','vector','BackgroundColor','none')
% %%
% print(Kdevol_fig,'iso_ana_num_comp_evol2','-dsvg','-r300')