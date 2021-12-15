%%
clc
clear all
close all

%% 3-DOF Simulation
global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tau1 tau2 tau3
L1 = 0.5; L2 = 0.5; L3 = 0.5;           % [m], link length
r1 = 0.1; r2 = 0.1; r3 = 0.1;           % [m], center of mass
m1 = 0.2; m2 = 0.2; m3 = 0.2;           % [kg], link mass
Iz1 = 0.05; Iz2 = 0.05; Iz3 = 0.05;     % [kgm^2], link inertia

g = 9.806;                              % [m/s^2], gravitational acceleration

% Initial pose
%    [rad]  [rad/s]
q1 = -pi/2; dq1 = 0;
q2 = pi/2;  dq2 = 0;
q3 = -pi/2; dq3 = 0;

%% Initial plot
FG = figure('Color', [1 1 1]);
AX = axes('parent', FG);
hold on;
grid on;
axis([-1.5 1.5 -1.5 1.5]);

x1 = L1*cos(q1);            % end coordinate x of first link
y1 = L1*sin(q1);            % end coordinate y of first link
Px1 = [0, x1];
Py1 = [0, y1];

x2 = L2*cos(q1+q2);         % end coordinate x of second link
y2 = L2*sin(q1+q2);         % end coordinate y of second link
Px2 = [x1, x1+x2];
Py2 = [y1, y1+y2];

x3 = L3*cos(q1+q2+q3);      % end coordinate x of third link
y3 = L3*sin(q1+q2+q3);      % end coordinate y of third link
Px3 = [x1+x2, x1+x2+x3];
Py3 = [y1+y2, y1+y2+y3];

p1 = plot(Px1, Py1, '-ob', 'Linewidth', 3);
p2 = plot(Px2, Py2, '-or', 'Linewidth', 3);
p3 = plot(Px3, Py3, '-og', 'Linewidth', 3);

xlabel('X-axis (m)', 'fontsize', 20);
ylabel('Y-axis (m)', 'fontsize', 20);
title('3-DOF Robot - Free fall', 'fontsize', 25);

%% ode45
% Simulation time
dt = 0.005; ft = 5;
n = 1;

% save as gif
filename = 'HW1_3_DOF_simulation.gif';

for cnt=0:dt:ft
    % Free fall simulation
    % The torque of each link is set to zero.
    tau1 = 0;               % Torque on the first link
    tau2 = 0;               % Torque on the second link
    tau3 = 0;               % Torque on the third link
    
    % Return the dynamics of each link
    [t,y] = ode45('three_link', [0 dt], [q1; dq1; q2; dq2; q3; dq3]);
    index = length(y);
    
    % Update dynamics of each link
    q1 = y(index, 1);
    dq1 = y(index, 2);
    q2 = y(index, 3);
    dq2 = y(index, 4);
    q3 = y(index, 5);
    dq3 = y(index, 6);
    
    % Update end point of first link
    x1 = L1*cos(q1);
    y1 = L1*sin(q1);
    Px1 = [0, x1];
    Py1 = [0, y1];
    
    % Update end point of second link
    x2 = L2*cos(q1+q2);
    y2 = L2*sin(q1+q2);
    Px2 = [x1, x1+x2];
    Py2 = [y1, y1+y2];
    
    % Update end point of third link
    x3 = L3*cos(q1+q2+q3);
    y3 = L3*sin(q1+q2+q3);
    Px3 = [x1+x2, x1+x2+x3];
    Py3 = [y1+y2, y1+y2+y3];
    
    % Print run time
    cmd = sprintf("Time: %2.2f", cnt);
    clc
    disp(cmd);
    
    % Set and plot current state
    n = n + 1;
    set(p1, 'XData', Px1, 'YData', Py1);
    set(p2, 'XData', Px2, 'YData', Py2);
    set(p3, 'XData', Px3, 'YData', Py3);
    drawnow
    
    % save as gif
%     frame = getframe(FG);
%     img = frame2im(frame);
%     [imind, cm] = rgb2ind(img, 256);
%     if cnt==0
%         imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1/60);
%     else
%         imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/60);
%     end
end
