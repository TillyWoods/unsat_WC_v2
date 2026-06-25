x = linspace(0,40,1000);
figure(1)
clf
plot(x,tanhAtoB(x,50,100,20,0.1))

[Qsi,Q0] = periodic_switching_forcing_smooth(dim);
t = linspace(0,120,1000);
figure(2)
clf
hold on
plot(t,Qsi(t))
%plot(t,Q0(t))

function Q = tanhAtoB(x,A,B,x0,eps)
    % tanh function that smoothly transitions from value Q=A to Q=B
    % over a window of size eps around x=x0.
    Q = ((B-A)/2)*tanh((x-x0)/eps) + (A+B)/2;
end