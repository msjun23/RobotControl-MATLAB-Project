%%
clc
clear all
close all

%% Dynamics parameter estimation
syms q dq ddq

global I Im m g r Fs Fv tau

I = 0.05;
Im = 0.05;
m = 0.5;
r = 0.2;
g = 9.806;
Fs = 0.1;
Fv = 0.1;

q = pi/4; dq = 0;

W1_int = [0 0 0 0];
theta = [0; 0; 0; 0];
P = eye(4);
u = 0;

data = [];
time = [];
n = 1;

%% Simulation
% Simulation time
dt = 0.005; ft = 5;

temp = zeros(4,4);
temp2 = zeros(4,1);
for cnt=0:dt:ft
    tau = sin(cnt) + cos(10*cnt);
    
    [t, y] = ode45('one_link', [0 dt], [q; dq]);
    index = length(y);
    
    q = y(index, 1);
    dq = y(index, 2);
    
    % W1 integration
    W1_int = W1_int + [0 -g*sin(q) -sign(dq) -dq]*dt;
    W2 = [dq 0 0 0];
    Y = W2 - W1_int;    % (1,4)
    u = u + tau*dt;     % (1x1)
    
    % Kalman Filter based parameter estimation algorithm
    P = P - P*Y.'*inv(eye(1) + Y*P*Y.')*Y*P;    % (4x4)
    K = P*Y.';                                  % (4x1)
    theta = theta + K*(u - Y*theta);            % (4x1)
    
    % Error minimization
%     temp = temp + Y.'*Y;
%     temp2 = temp2 + Y.'*u;
%     theta = inv(temp)*temp2;
    
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

plot(time(:,1), data(:,1), 'r', 'linew', 2);
grid on;
hold on;
plot(time(:,1), data(:,2), 'g', 'linew', 2);
plot(time(:,1), data(:,3), 'b', 'linew', 2);
plot(time(:,1), data(:,4), 'k', 'linew', 2);

xlabel('time', 'fontsize', 15);
ylabel('parameter', 'fontsize', 15);
legend({'I+Im=0.1', 'mr=0.1', 'Fs=0.1', 'Fv=0.1'}, 'Location', 'southeast');
