%% Define
function G = GetGravity(q)
global g L;

G = g/L*sin(q);
