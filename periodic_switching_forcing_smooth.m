function [Qsi,Q0] = periodic_switching_forcing_smooth(dim)
    % Periodic switching forcing. This consists of 3 phases:
    % - growth phase
    % - freezing phase
    % - melt phase
    % with sudden switches between each. The output is shotrwave radiation
    % Qsi and longwave radition and turbulent heat fluxes Q0, both as
    % function handles of t.

    dim.td = 1;

    % Growth phase
    Qsi1 = 100;
    Q01 = 50;
    gap1 = 20;      % duration of section
    
    % Freeze phase
    Qsi2 = 0;
    Q02 = -50;
    gap2 = 20;      % duration of section
    
    % Melt phase
    Qsi3 = 0;
    Q03 = 100;
    gap3 = 20;      % duration of section
    
    period = gap1 + gap2 + gap3;    % total period

    eps = 0.1;      % time over which to smoothly transition between phases

    tanhQsi0 = @(t)tanhAtoB(t,Qsi3,Qsi1,0,eps*dim.td);
    tanhQsi1 = @(t)tanhAtoB(t,Qsi1,Qsi2,gap1*dim.td,eps*dim.td);
    tanhQsi2 = @(t)tanhAtoB(t,Qsi2,Qsi3,(gap1+gap2)*dim.td,eps*dim.td);
    tanhQsi3 = @(t)tanhAtoB(t,Qsi3,Qsi1,period*dim.td,eps*dim.td);

    tanhQ00 = @(t)tanhAtoB(t,Q03,Q01,0,eps*dim.td);
    tanhQ01 = @(t)tanhAtoB(t,Q01,Q02,gap1*dim.td,eps*dim.td);
    tanhQ02 = @(t)tanhAtoB(t,Q02,Q03,(gap1+gap2)*dim.td,eps*dim.td);
    tanhQ03 = @(t)tanhAtoB(t,Q03,Q01,period*dim.td,eps*dim.td);

    
    % Periodic functions for Qsi and Q0
    Qsi = @(t) tanhQsi0(mod(t,period*dim.td)).*and(mod(t,period*dim.td)<=0.5*gap1*dim.td,t>0.5*gap1*dim.td)... 
        + tanhQsi1(mod(t,period*dim.td)).*and(mod(t,period*dim.td)>0.5*gap1*dim.td,mod(t,period*dim.td)<=(gap1+0.5*gap2)*dim.td)...
        + tanhQsi2(mod(t,period*dim.td)).*and(mod(t,period*dim.td)>(gap1+0.5*gap2)*dim.td,mod(t,period*dim.td)<=(gap1+gap2+0.5*gap3)*dim.td)...
        + tanhQsi3(mod(t,period*dim.td)).*(mod(t,period*dim.td)>(gap1+gap2+0.5*gap3)*dim.td)...
        + Qsi1*(t<=0.5*gap1*dim.td);

    Q0 = @(t) tanhQ00(mod(t,period*dim.td)).*and(mod(t,period*dim.td)<=0.5*gap1*dim.td,t>0.5*gap1*dim.td)...
        + tanhQ01(mod(t,period*dim.td)).*and(mod(t,period*dim.td)>0.5*gap1*dim.td,mod(t,period*dim.td)<=(gap1+0.5*gap2)*dim.td)...
        + tanhQ02(mod(t,period*dim.td)).*and(mod(t,period*dim.td)>(gap1+0.5*gap2)*dim.td,mod(t,period*dim.td)<=(gap1+gap2+0.5*gap3)*dim.td)...
        + tanhQ03(mod(t,period*dim.td)).*(mod(t,period*dim.td)>(gap1+gap2+0.5*gap3)*dim.td)...
        + Q01*(t<=0.5*gap1*dim.td);


    function Q = tanhAtoB(x,A,B,x0,eps)
        % tanh function that smoothly transitions from value Q=A to Q=B
        % over a window of size eps around x=x0.
        Q = ((B-A)/2)*tanh((x-x0)/eps) + (A+B)/2;
    end

end