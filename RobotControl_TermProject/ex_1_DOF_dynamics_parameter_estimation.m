%% Dynamics parameter estimation
clc
clear all
close all

syms q dq ddq

global I Im m g r Fs Fv tau

I = 0.05;
Im = 0.05;
m = 0.5;
r = 0.2;
g = 9.806;
Fs = 0.1;
Fv = 0.1;

dt = 0.002; ft = 5;
q = pi/4; dq = 0;

W1_int = [0 0 0 0];
theta = [0; 0; 0; 0];
P = eye(4);
u = 0;

data = [];
time = [];
n = 1;

for cnt=0:dt:ft
    tau = sin(cnt) + cos(10*cnt);
    
    [t, y] = ode45('one_link', [0 dt], [q; dq]);
    index = length(y);
    
    q = y(index, 1);
    dq = y(index, 2);
    
    % Kalman Filter based parameter estimation algorithm
    W1_int = W1_int + [0 -g*sin(q) -sign(dq) -dq]*dt;
    W2 = [dq 0 0 0];
    Y = W2 - W1_int;
    u = u + tau*dt;
    P = P - P*Y.'*inv(eye(1) + Y*P*Y.')*Y*P;
    K = P*Y.';
    theta = theta + K*(u - Y*theta);
    
    time(n,:) = cnt;
    data(n,:) = theta;
    n = n + 1;
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
legend({'I + Im', 'm * r', 'Fs', 'Fv'}, 'Location', 'southeast');
