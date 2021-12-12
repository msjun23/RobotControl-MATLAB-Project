%% Define
function X = GetKinematics(q1, q2)
global L1 L2

x = L1*cos(q1) + L2*cos(q1+q2);
y = L1*sin(q1) + L2*sin(q1+q2);

X = [x; y];
