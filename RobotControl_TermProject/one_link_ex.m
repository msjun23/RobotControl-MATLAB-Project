%% Define
function dydt = one_link_ex(t, y)
global m L g tau;

dydt = [y(2); -g/L*sin(y(1)) + tau/m/L/L];
