%% Define
function dydt = one_link_ex(t, y)
global m I g u;

dydt = [y(2); -g/I*sin(y(1)) + u/m/I/I];
