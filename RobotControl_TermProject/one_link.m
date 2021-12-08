%% Define
function dydt = one_link(t, y)
global I Im m g r Fs Fv tau

ddy = (tau - m*g*r*sin(y(1)) - Fs*sign(y(2)) - Fv*y(2))/(I + Im);

dydt = [y(2); ddy];
