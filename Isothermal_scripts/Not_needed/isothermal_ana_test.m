% Isothermal water transport analytical test
% 23rd May 2024

t = linspace(0,1,11);
Z = linspace(0,1);

% % Parameters
% dim = dimensional_parameters();
% dim.Kd = 1;
% dim.Ld = 100;
% scale = scalings(dim);
% pp = dimensionless_parameters(dim,scale);
% pp = isothermal_parameters(pp);

phi0 = 0.2;
Z0 = 0.5;
q0 = 0.1;
K = 100;

Kdmax = 3/Z0;

Kd = 2;

[S,q,qd,phi,Zb_func,Zt_func,Zw,ts]=isothermal_drainage_analytic(t,Z,K,Kd,phi0,Z0,q0);

figure(3)
clf
tiledlayout(2,5)
for i=1:10
    nexttile
    box on
    hold on
    p1 = plot(phi,Z,'DisplayName','$\phi$');
    p2 = plot(S(:,i),Z,'DisplayName','$S$');
    plot(10*q(:,i),Z,'DisplayName','$10q$');
    set(gca,'YDir','reverse')
    if or(i==1,i==6)
        ylabel('Z')
    end
    % legend('show','interpreter','latex')

end

%% Saturation time dependence of Kd
Zw = @(t,Kd) (3./Kd).*log(1 + (1/3).*Kd.*K^(1/3).*q0^(2/3).*t);
% Point at which saturation first occurs
ts = @(Kd)(3./(Kd.*K^(1/3).*q0^(2/3))).*((phi0.*(K/q0)^(1/3)).^(Kd./(3/Z0-Kd))-1);
Zs = @(Kd) Zw(ts(Kd),Kd);

Kd = linspace(0,Kdmax);
figure(2)
clf
box on
hold on
plot(Kd,ts(Kd),'DisplayName','$t_\mathrm{s}$')
plot(Kd,Zs(Kd),'DisplayName','$Z_\mathrm{s}$')
xlabel('$\mathcal{K}_\mathrm{d}$')
legend('show','interpreter','latex','location','northwest')
ylim([0,20])

%% Wetting front drainage and no drainage comparison
Zw_nd = @(t) q0^(2/3)*K^(1/3)*t;

Kd_vals = 0:10;
t = linspace(0,5);

figure(2)
clf
box on
hold on
for i = 1:length(Kd_vals)
    Kd = Kd_vals(i);
    if Kd==0
        plot(t,Zw_nd(t),'DisplayName','$\mathcal{K}_\mathrm{d} = 0$')
    else
        plot(t,Zw(t,Kd),'DisplayName',['$\mathcal{K}_\mathrm{d} = $',num2str(Kd)])
    end
end
set(gca,'YDir','reverse')
xlabel('$t$')
ylabel('$Z_\mathrm{w}$')
legend('show','interpreter','latex','location','northwest')

%% Saturation front drainage and ni drainage comparison
Kd_vals = 0:6;
t = linspace(0,10,10000);
Z = linspace(0,10,1000);

cmapS = parula(length(Kd_vals));

phi0 = 0.2;
Z0 = 0.5;
q0 = 0.1;
K = 100;

figure(10)
clf
box on
hold on
for i = 1:length(Kd_vals)
    Kd = Kd_vals(i);
    [S,q,phi,Zb_func,Zt_func,Zw,ts]=isothermal_drainage_analytic(t,Z,K,Kd,phi0,Z0,q0);
    sat_surf = find(S(1,:)==1,1,'first');
    sat_surf(isempty(sat_surf)) = length(t);
    tsat = t(sat_surf);
    plot(t(t<tsat)-ts,Zt_func(t(t<tsat))-Zw(ts),'DisplayName',['$\mathcal{K}_\mathrm{d} = $',num2str(Kd)],'color',cmapS(i,:))
    plot(t(t<tsat)-ts,Zb_func(t(t<tsat))-Zw(ts),'DisplayName',['$\mathcal{K}_\mathrm{d} = $',num2str(Kd)],'color',cmapS(i,:))
end
set(gca,'YDir','reverse')
xlabel('$t-t_\mathrm{s}$')
ylabel('Saturation fronts')
xlim([0,1])
ylim([-1,1])
% legend('show','interpreter','latex','location','northwest')

%% Colour plots of solutions
t = linspace(0,5,500);
Z = linspace(0,3,100);

phi0 = 0.2;
Z0 = 0.5;
q0 = 0.1;
K = 100;

Kdmax = 3/Z0;

Kd = 2;

dratio = 3;

[S,q,qd,phi,Zb_func,Zt_func,Zw,ts]=isothermal_drainage_analytic(t,Z,K,Kd,phi0,Z0,q0);

sat_surf = find(S(1,:)==1,1,'first');
sat_surf(isempty(sat_surf)) = length(t);
tsat = t(sat_surf);

Zend = Z(end);
tend = t(sat_surf);

h = 6;

figure(14)
clf
box on
hold on
set(gcf,'units','centimeters','paperpositionmode','auto','position',[10 10 h*dratio*tend/Zend+1 h],'render','painters')
[tgrid,Zgrid] = meshgrid(t,Z);
surf(tgrid,Zgrid,S,'LineStyle','none')
view(0,90);
line(t(t<ts),Zw(t(t<ts)),2*t(t<ts).^0,'linewidth',2,'color','red','linestyle','-')
if ~isinf(ts)
line(t(t>ts),Zt_func(t(t>ts)),2*t(t>ts).^0,'linewidth',2,'color','k','linestyle','-')
line(t(t>ts),Zb_func(t(t>ts)),2*t(t>ts).^0,'linewidth',2,'color',['k'],'linestyle','-')
end
xlim([min(t),tend])
ylim([0,Zend])
xlabel('$t$','interpreter','latex')
ylabel('$Z$','interpreter','latex')
title(['$\mathcal{K}_\mathrm{d}=',num2str(Kd),'$'],'interpreter','latex')
set(gca,'YDir','reverse')
c1 = colorbar;
caxis([0,1])
ylabel(c1,'$S$','interpreter','latex')
% daspect([1,dratio,1])
% text(0.2,0.25,'$S=0$','units','normalized','color',[1 1 1])
% text(0.7,0.6,'$S=1$','units','normalized')
% fig_name = 'isothermal_Sana';
% print(gcf,'-depsc2',fig_name,'-loose');

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

figure(12)
clf
box on
hold on
set(gcf,'units','centimeters','paperpositionmode','auto','position',[10 10 18 5]) %,'render','painters')
% tiledlayout(1,3,'padding','compact')

h = 0.9;

y0 = 0.07;
gaps = 0.02;
widths = [0.1,0.26,0.35+0.09]*1.08;
widths = [0.1,0.26,0.26+0.073]*1.25;
widths = [0.1,0.26*4/2.2,(0.26+0.07)*4/2.2]*1.25*0.48*1.25;
x0 = [2*gaps,3*gaps+widths(1),4*gaps+widths(2)+widths(1)]+gaps;
w = [widths]; %0.2,0.3,0.4];

clear ax;

for i=1:length(Kd_vals)
    
    Kd = Kd_vals(i);
[S,q,phi,Zb_func,Zt_func,Zw,ts]=isothermal_drainage_analytic(t,Z,K,Kd,phi0,Z0,q0);

% w = h*Zend/ts/dratio;
ax(i) = axes('position',[x0(i),y0,w(i),h],'units','normalized');
% nexttile
% x0 = x0+w;

sat_surf = find(S(1,:)==1,1,'first');
sat_surf(isempty(sat_surf)) = length(t);

Zend = 2; %Z(end);
tend = min(t(sat_surf),4); %2.2);

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
xlabel('$t$','interpreter','latex')
if i==1
ylabel('$Z$','interpreter','latex')
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

end
% text(0.2,0.25,'$S=0$','units','normalized','color',[1 1 1])
% text(0.7,0.6,'$S=1$','units','normalized')
% fig_name = 'isothermal_Sana';
% print(gcf,'-depsc2',fig_name,'-loose');
