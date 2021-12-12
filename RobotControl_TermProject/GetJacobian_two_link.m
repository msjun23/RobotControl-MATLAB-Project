%% Define
function J = GetJacobian_two_link(q1, q2)
global L1 L2;

J = [-L1*sin(q1)-L2*sin(q1+q2) -L2*sin(q1+q2); 
     L1*cos(q1)+L2*cos(q1+q2)  L2*cos(q1+q2)];
