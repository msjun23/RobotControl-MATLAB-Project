%%
clc
clear all
close all

%% Dynamics parameter estimation
syms q1 q2 dq1 dq2 ddq1 ddq2

global I1 I2 Im1 Im2 L1 L2 m1 m2 r1 r2 g Fs1 Fs2 Fv1 Fv2 tau1 tau2

I1 = 0.05;	I2 = 0.05;
Im1 = 0.05;	Im2 = 0.05;
L1 = 0.5;	L2 = 0.5;
m1 = 0.2;	m2 = 0.2;
r1 = 0.1;	r2 = 0.1;
g = 9.806;
Fs1 = 0.1;	Fs2 = 0.1;
Fv1 = 0.1;	Fv2 = 0.1;

q1 = pi/4;  dq1 = 0;
q2 = pi/4; dq2=0;

W1_int = [0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0];
theta = [0; 0; 0; 0; 0;
         0; 0; 0; 0; 0];
P = eye(10);
u = [0; 0];

data = [];
time = [];
n = 1;

%% Simulation
% Simulation time
dt = 0.005; ft = 5;

temp = zeros(10,10);
temp2 = zeros(10,1);
for cnt=0:dt:ft
    tau1 = sin(cnt) + cos(10*cnt);
    tau2 = sin(cnt) + cos(10*cnt);
    
    [t, y] = ode45('HW2_two_link', [0 dt], [q1; dq1; q2; dq2]);
    index = length(y);
    
    q1 = y(index, 1);
    dq1 = y(index, 2);
    q2 = y(index, 3);
    dq2 = y(index, 4);
    
    % W1 integration
    W1_int = W1_int + [0 0 0                      g*cos(q1) g*cos(q1+q2) 0 -sign(dq1) 0          -dq1 0; 
                       0 0 -dq1*(dq1+dq2)*sin(q2) 0         g*cos(q1+q2) 0 0          -sign(dq2) 0    -dq2]*dt;
    W2 = [dq1 dq1+dq2 (2*dq1+dq2)*cos(q2) 0 0 0   0 0 0 0;
          0   dq1+dq2 dq1*cos(q2)         0 0 dq2 0 0 0 0];
    Y = W2 - W1_int;            % (2x10)
    u = u + [tau1; tau2]*dt;    % (2x1)
    
    % Kalman Filter based parameter estimation algorithm
    P = P - P*Y.'*inv(eye(2) + Y*P*Y.')*Y*P;    % (10x10)
    K = P*Y.';                                  % (10x1)
    theta = theta + K*(u - Y*theta);            % (10x1)
    
    % Error minimization
    temp = temp + Y.'*Y;
    temp2 = temp2 + Y.'*u;
    theta = inv(temp)*temp2;
    
    time(n,:) = cnt;
    data(n,:) = theta;
    n = n + 1;
    
    cmd = sprintf("Time: %2.2f", cnt);
    clc
    disp(cmd);
end

%% Plot
FG = figure('Color', [1 1 1]);
AX = axes('parent', FG);

plot(time(:,1), data(:,1), 'k', 'linew', 2);
grid on;
hold on;
plot(time(:,1), data(:,2), 'r', 'linew', 2);
plot(time(:,1), data(:,3), 'g', 'linew', 2);
plot(time(:,1), data(:,4), 'b', 'linew', 2);
plot(time(:,1), data(:,5), 'y', 'linew', 2);
plot(time(:,1), data(:,6), '--k', 'linew', 2);
plot(time(:,1), data(:,7), 'c', 'linew', 2);
plot(time(:,1), data(:,8), '--c', 'linew', 2);
plot(time(:,1), data(:,9), 'm', 'linew', 2);
plot(time(:,1), data(:,10), '--m', 'linew', 2);

xlabel('time', 'fontsize', 15);
ylabel('parameter', 'fontsize', 15);
legend({'I1+m2L1^2+Im1=0.15', 'I2=0.05', 'm2r2L1=0.01', 'm1r1+m2L1=0.12', 'm2r2=0.02', ...
        'Im2=0.05', 'Fs1=0.1', 'Fs2=0.1', 'Fv1=0.1', 'Fv2=0.1'}, ...
        'Location', 'best');
