function [Qsi,Q0] = periodic_switching_forcing(dim)
    % Periodic switching forcing. This consists of 3 phases:
    % - growth phase
    % - freezing phase
    % - melt phase
    % with sudden switches between each. The output is shotrwave radiation
    % Qsi and longwave radition and turbulent heat fluxes Q0, both as
    % function handles of t.

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
    
    % Periodic functions for Qsi and Q0
    Qsi = @(t) Qsi1*(mod(t,period*dim.td)<gap1*dim.td)...
        +Qsi2*and(mod(t,period*dim.td)>=gap1*dim.td,mod(t,period*dim.td)<=(gap1+gap2)*dim.td)...
        +Qsi3*(mod(t,period*dim.td)>(gap1+gap2)*dim.td);

    Q0 = @(t) Q01*(mod(t,period*dim.td)<gap1*dim.td)...
        +Q02*and(mod(t,period*dim.td)>=gap1*dim.td,mod(t,period*dim.td)<=(gap1+gap2)*dim.td)...
        +Q03*(mod(t,period*dim.td)>(gap1+gap2)*dim.td);

end