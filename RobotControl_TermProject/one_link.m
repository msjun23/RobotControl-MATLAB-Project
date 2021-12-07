%% Define
function dydt = one_link(t, y)
% global l lm m g r Fs Fv tau
% 
% ddy = (tau - m*g*r*sin(y(1)) - Fs*sign(y(2)) - Fv*y(2))/(l + lm);
% 
% dydt = [y(2); ddy];

global m l g u;

dydt = [y(2); -g/l*sin(y(1)) + u/m/l/l];
