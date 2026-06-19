function [S,q,qd,phi,Zb_func,Zt_func,Zw,ts]=isothermal_drainage_analytic(t,Z,K,Kd,phi0,Z0,q0)
% Isothermal water transport WITH lateral drainage - semi-analytical solution
% 22nd May 2024 (modified from no drainage case - 12th April 2021)
%
% Neglecting capillary pressure.
% A constant flux q0 is prescribed at the surface, and the poroisty profile
% is prescribed as phi(Z) = phi0*exp(-Z/Z0).
% K is the dimensionles permeability, Kd is the dimensionless lateral
% drainage strength
%
% Zw is the wetting front (before saturation is reached).
% Zt and Zb are the top and bottom saturation fronts.
% ts is the time at which saturation first occurs.

tend=t(end);

% Prescribed porosity
phi = @(Z) phi0*exp(-Z/Z0);

if Kd == 0  % Solution if there is no lateral drainage
    % Wetting front (before saturation is reached)
    Zw = @(t) q0^(2/3)*K^(1/3)*t;

    % Saturation and flux in unsaturated region
    Su = @(Z) (1./phi(Z)).*(q0/K)^(1/3);
    qu = @(Z) q0;
    
    % Point at which saturation first occurs
    ts = -(Z0/(3*q0^(2/3)*K^(1/3)))*log(q0/(K*phi0^3));
    Zs = Zw(ts);
    
    % Saturation fronts (top Zt and bottom Zb of saturated region)
    % First solve for early times using asymptotics in order to avoid the
    % intial signularity in the numerics
    epsilon = 1e-4;     % small time
    Wt=(1/4)*(sqrt(33)-3)*K^(1/3)*q0^(2/3);     % corresponding O(e) change in Zt
    Wb=K^(1/3)*q0^(2/3);                        % corresponding O(e) chnge in Zb
    qs = @(Zt,Zb,Z) ((3.*K*(phi0^3))./Z0).*(Zt-Zb)./(exp(3.*Zt./Z0)-exp(3.*Zb./Z0));    % (constant) flux in saturated region)
    func = @(t,y) [(qs(y(1),y(2))-q0)./(phi(y(1)).*(1-(q0./(K.*(phi(y(1))^3)))^(1/3))); qs(y(1),y(2))./phi(y(2))];  % RHS of ODE
    y0 = [Zs-epsilon*Wt,Zs+epsilon*Wb];         % inital condition (both top and bottom fronts are at the wetting front)
    [tsat,y]=ode15s(func,[ts+epsilon,tend],y0);
    Zt=y(:,1);
    Zb=y(:,2);
    Zt_func = @(t) interp1(tsat,Zt,t);      % function to calculate Zt for any time t>ts
    Zb_func = @(t) interp1(tsat,Zb,t);      % function to calculate Zb for any time t>ts

else    % Solution WITH lateral drainage
    % Wetting front (before saturation is reached)
    % Zw = @(t) (3/(Kd))*log(1 + (1/3)*Kd*K^(1/3)*q0^(2/3)*t);  % WRONG!!!
    Zw = @(t) (3/(2*Kd))*log(1 + (2/3)*Kd*K^(1/3)*q0^(2/3)*t);
    
    % Saturation and flux in unsaturated region
    Su = @(Z) (1./phi(Z)).*(q0/K)^(1/3).*exp(-Kd*Z/3);
    qu = @(Z) q0*exp(-Kd*Z);
    
    if Kd<3/Z0
    % Point at which saturation first occurs
    % ts = (3/(Kd*K^(1/3)*q0^(2/3)))*((phi0*(K/q0)^(1/3))^(Kd/(3/Z0-Kd))-1); % WRONG
    ts = (3/(2*Kd*K^(1/3)*q0^(2/3)))*((phi0*(K/q0)^(1/3))^(2*Kd/(3/Z0-Kd))-1);
    Zs = Zw(ts);
    
    % Saturation fronts (top Zt and bottom Zb of saturated region)
    % First solve for early times using asymptotics in order to avoid the
    % intial singularity in the numerics
    epsilon = 1e-4;     % small time
    Wt=1; %(1/4)*(sqrt(33)-3)*K^(1/3)*q0^(2/3);     % corresponding O(e) change in Zt
    Wb=1; %K^(1/3)*q0^(2/3);                        % corresponding O(e) chnge in Zb
    % qs0 = @(Zt,Zb) (3/Z0 - Kd)*K*(phi0^3)*(Zb-Zt)/(exp((3/Z0-Kd)*Zb)-exp((3/Z0-Kd)*Zt));    % constant used i flux in saturated region
    qs0 = @(Zt,Zb) K*(phi0^3)*(3/Z0-Kd)*(Zb-Zt)./(exp(3*Zb/Z0)-exp(3*Zt/Z0));
    % qs = @(Zt,Zb,Z) qs0(Zt,Zb)*exp(-Kd*Z);    % flux in saturated region
    qs = @(Zt,Zb,Z) qs0(Zt,Zb) + (Z0/3)*Kd*K*phi(Z).^3;    % flux in saturated region
    % func = @(t,y) [(qs0(y(1),y(2))-q0).*exp(-Kd*y(1))./(phi(y(1))-(q0/K)^(1/3)*exp(-Kd*y(1)/3));...
    %     qs0(y(1),y(2)).*exp(-Kd*y(2))./phi(y(2))];  % RHS of ODE
    func = @(t,y) [(qs(y(1),y(2),y(1))-qu(y(1)))./(phi(y(1))-(q0/K)^(1/3).*exp(-Kd*y(1)/3));...
        qs(y(1),y(2),y(2))./phi(y(2))];
    y0 = [Zs-epsilon*Wt,Zs+epsilon*Wb];         % inital condition (both top and bottom fronts are at the wetting front)
    [tsat,y]=ode15s(func,[ts+epsilon,tend],y0);
    Zt=y(:,1);
    Zb=y(:,2);
    Zt_func = @(t) interp1(tsat,Zt,t);      % function to calculate Zt for any time t>ts
    Zb_func = @(t) interp1(tsat,Zb,t);      % function to calculate Zb for any time t>ts
    else
        ts = Inf;
        Zs = NaN;
        Zt_func = @(t) NaN;
        Zb_func = @(t) NaN;
    end
end

S=zeros(size(Z));
q=zeros(size(Z));

Sstore=zeros(length(Z),length(t));
qstore=zeros(length(Z),length(t));
qdstore=zeros(length(Z),length(t));

for i=[1:length(t)]
    if t(i)<=ts
        above = (Z<=Zw(t(i)));
        below = (Z>Zw(t(i)));
        S(above) = Su(Z(above));    % saturation above wetting front (unsaturated)
        S(below) = 0;                             % no water below wetting front
        q(above) = qu(Z(above));               % water flux above wetting front (unsaturated)
        q(below) = 0;                             % no water flux below
    end
    if t(i)>ts
        above = (Z<Zt_func(t(i)));
        below = (Z>Zb_func(t(i)));
        sat = logical(1-above-below);
        S(above) = Su(Z(above));     % saturation as before in unsaturated region
        S(sat) = 1;                               % S=1 in saturated region
        S(below) = 0;                             % no water below saturated region
        q(above) = qu(Z(above));                            % unsaturated flux above
        q(sat) = qs(Zt_func(t(i)),Zb_func(t(i)),Z(sat));                 % saturated flux in saturated region
        q(below) = 0;                             % no water flux below
    end
    qd = Kd*K*phi(Z).^3.*S.^3;   % lateral drainage
    qdstore(:,i)=qd;
    Sstore(:,i)=S;
    qstore(:,i)=q;
end

S=Sstore;
q=qstore;
phi=phi(Z);
qd = qdstore;
end