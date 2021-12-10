%% Define
function G = GetGravity(th)
global I g L;

G = g/L*sin(th);
