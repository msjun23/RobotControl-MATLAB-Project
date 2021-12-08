%% 2-DOF Simulation
clear all
close all
global Iz1 Iz2 L1 L2 g m1 m2 r1 r2 tau1 tau2

L1 = 0.5; L2 = 0.5; L3 = 0.5;
r1 = 0.1; r2 = 0.1; r3 = 0.1;
m1 = 0.2; m2 = 0.2; m3 = 0.2;
Iz1 = 0.05; Iz2 = 0.05; Iz3 = 0.05;

g = 9.806;

dt = 0.005; ft = 5;
q1 = -pi/2; dq1 = 0;
q2 = pi/2; dq2 = 0;

data = [];
n = 1;

FG = figure('Color', [1 1 1]);
AX = axes('parent', FG);
hold on;
grid on;
axis([-1.5 1.5 -1.5 1.5]);

x1 = L1*cos(q1);
y1 = L1*sin(q1);
Px1 = [0, x1];
Py1 = [0, y1];

x2 = L2*cos(q1+q2);
y2 = L2*sin(q1+q2);
Px2 = [x1, x1+x2];
Py2 = [y1, y1+y2];

p1 = plot(Px1, Py1, '-ob', 'Linewidth', 3);
p2 = plot(Px2, Py2, '-or', 'Linewidth', 3);

%% ode45
for cnt=0:dt:ft
    tau1 = 0;
    tau2 = 0;
    
    [t,y] = ode45('two_link', [0 dt], [q1; dq1; q2; dq2]);
    
    index = length(y);
    
    q1 = y(index, 1);
    dq1 = y(index, 2);
    q2 = y(index, 3);
    dq2 = y(index, 4);
    
    x1 = L1*cos(q1);
    y1 = L1*sin(q1);
    Px1 = [0, x1];
    Py1 = [0, y1];
    
    x2 = L2*cos(q1+q2);
    y2 = L2*sin(q1+q2);
    Px2 = [x1, x1+x2];
    Py2 = [y1, y1+y2];
    
    n = n + 1;
    cmd = sprintf("Time: %2.2f", cnt);
    clc
    disp(cmd);
    
    if rem(n, 10) == 0
        set(p1, 'XData', Px1, 'YData', Py1);
        set(p2, 'XData', Px2, 'YData', Py2);
        drawnow
    end
end
