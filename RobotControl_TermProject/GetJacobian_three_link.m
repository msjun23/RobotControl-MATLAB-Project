%% Define
function J = GetJacobian_three_link(q1, q2, q3)
global L1 L2 L3;

J = [-L1*sin(q1)-L2*sin(q1+q2)-L3*sin(q1+q2+q3) -L2*sin(q1+q2)-L3*sin(q1+q2+q3) -L3*sin(q1+q2+q3);
     L1*cos(q1)+L2*cos(q1+q2)+L3*cos(q1+q2+q3)  L2*cos(q1+q2)+L3*cos(q1+q2+q3)  +L3*cos(q1+q2+q3)];
