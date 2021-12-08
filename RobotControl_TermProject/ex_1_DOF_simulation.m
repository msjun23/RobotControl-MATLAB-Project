%% 1-DOF Simulation
global m l g u;

dt = 0.005; ft = 5;
q = pi/4; dq = 0;
data = [];

m = 1; l = 1; g = 9.8148; n =1; u = 0;

FG = figure('Color', [1 1 1]);
AX = axes('parent', FG);
hold on;
grid on;
axis([-1.5 1.5 -1.5 1.5]);

Px = [0, 1];
Py = [0, 0];

p = plot(Px, Py, '-ob', 'Linewidth', 3);

%% ODE45
for cnt=0:dt:ft
    [t,y] = ode45('one_link', [0 dt], [q; dq]);
    
    index = length(y);
    
    q = y(index, 1);
    dq = y(index, 2);
    
    x = l*sin(q);
    y = -l*cos(q);
    Px = [0 x];
    Py = [0 y];
    
    data(n,1) = cnt;
    data(n,2) = q;
    data(n,3) = dq;
    n = n + 1;
    
    cmd = sprintf("Time: %2.2f", cnt);
    clc
    disp(cmd);
    
    if rem(n, 10) == 0
        set(p, 'XData', Px, 'YData', Py);
        drawnow
    end
end

%%
FG2 = figure('Color', [1 1 1]);
AX2 = axes('parent', FG2);
grid on;

plot(data(:,1), data(:,2), 'r');
hold on;
plot(data(:,1), data(:,3), 'b');
