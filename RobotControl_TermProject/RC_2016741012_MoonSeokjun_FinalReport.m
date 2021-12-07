%% Clear
clc
clear all
close all

%% DH Parameter, Tranformation matrix

syms L1 L2 m1 m2 Ic1 Ic2 Im1 Im2 r1 r2 Iz1 Iz2
syms th1 th2 dth1 dth2 ddth1 ddth2

I1xx = 0;
I1yy = Iz1;
I1zz = Iz1;

I2xx = 0;
I2yy = Iz2;
I2zz = Iz2;

d1 = 0;     d2 = 0;
a1 = L1;    a2 = L2;
al1 = 0;    al2 = 0;

DH = [th1, d1, a1, al1 ;
      th2, d2, a2, al2];

T01 = HT(th1, d1, a1, al1);
T12 = HT(th2, d2, a2, al2);
T02 = T01 * T12;

%% Differential matrix(rotation)
Qr = [0 -1 0 0;
      1  0 0 0;
      0  0 0 0;
      0  0 0 0];

Q1 = Qr;
Q2 = Qr;

dT01t1 = Qr * T01;
% dT01t2 = diff(T01);       <- same result

%% U matrix
U11 = Q1 * T01;             % i=1, j=1
U12 = zeros(4, 4);          % i=1, j=2

U21 = Q1 * T02;             % i=2, j=1
U22 = T01 * Q2 * T12;       % i=2, j=2

%% Pseudo-inverse
J1(1, 1) = (-I1xx + I1yy + I1zz) / 2;
J1(2, 2) = ( I1xx - I1yy + I1zz) / 2;
J1(3, 3) = ( I1xx + I1yy - I1zz) / 2;
J1(1, 4) = -m1 * (L1 - r1);
J1(4, 1) = J1(1, 4);
J1(4, 4) = m1;

J2(1, 1) = (-I2xx + I2yy + I2zz) / 2;
J2(2, 2) = ( I2xx - I2yy + I2zz) / 2;
J2(3, 3) = ( I2xx + I2yy - I2zz) / 2;
J2(1, 4) = -m2 * (L2 - r2);
J2(4, 1) = J2(1, 4);
J2(4, 4) = m2;

%% Inertia matrix
D11 = trace(U11 * J1 * U11.')...    % i=1, k=1, j=1
    + trace(U21 * J2 * U21.');      % i=1, k=1, j=2
M11 = simplify(D11);

D12 = trace(U22 * J2 * U21.');      % i=1, k=2, j=2
M12 = simplify(D12);

M21 = M12;

D22 = trace(U22 * J2 * U22.');      % i=2, k=2, j=2
M22 = simplify(D22);

M = [M11 M12; M21 M22];

%% Inertia matrix -> Define using summation function
n = 2;
for i = 1:n
    for k = 1:n
        nM(i,k) = Inertia(i, k, n, U11, U12, U21, U22, J1, J2);
    end
end

%% Coriolis & Centrifugal matrix
U111 = Q1 * Q1 * T01;       % i=1, j=1, k=1
U112 = zeros(4, 4);         % i=1, j=1, k=2
U121 = zeros(4, 4);         % i=1, j=2, k=1
U122 = zeros(4, 4);         % i=1, j=2, k=2
U211 = Q1 * Q1 * T02;       % i=2, j=1, k=1
U212 = Q1 * T01 * Q2 * T12; % i=2, j=1, k=2
U221 = Q2 * T01 * Q2 * T12; % i=2, j=2, k=1
U222 = T01 * Q2 * Q2 * T12; % i=2, j=2, k=2

%%  Coriolis & Centrifugal matrix -> Define using function
n = 2;

for i = 1:n
    for j = 1:n
        for k = 1:n
            cmd = sprintf("nU%d%d%d = dUdq(i, j, k, T01, T12, T02, Q1, Q2);", i, j, k);
            eval(cmd);
        end
    end
end

%% h term
h111 = trace(U111 * J1 * U11.')...      % i=1, k=1, m=1, j=1
     + trace(U211 * J2 * U21.');        % i=1, k=1, m=1, j=2    dth1*dth1
h112 = trace(U212 * J2 * U21.');        % i=1, k=1, m=2, j=2    dth1*dth2
h121 = trace(U221 * J2 * U21.');        % i=1, k=2, m=1, j=2    dth2*dth1
h122 = trace(U222 * J2 * U21.');        % i=1, k=2, m=2, j=2    dth2*dth2

h1 = (dth1^2)*h111 + (dth1*dth2)*(h112 + h121) + (dth2^2)*h122;

h211 = trace(U211 * J2 * U22.');        % i=2, k=1, m=1, j=2
h212 = trace(U212 * J2 * U22.');        % i=2, k=1, m=2, j=2
h221 = trace(U221 * J2 * U22.');        % i=2, k=2, m=1, j=2
h222 = trace(U222 * J2 * U22.');        % i=2, k=2, m=2, j=2

h2 = (dth2^2)*h211 + (dth1*dth2)*(h212 + h221) + (dth2^2)*h222;

h = simplify([h1; h2]);

%%
L2 = L1;

r1 = L1 / 2;
r2 = L2 / 2;
Iz1 = 1/3 * m1 * L1^2;
Iz2 = 1/3 * m2 * L2^2;

sh1 = simplify(eval(h1));
sh2 = simplify(eval(h2));

%% Gravity term
n = 2;
syms r1 r2 Iz1 Iz2
syms g
r11 = [-(L1-r1); 0; 0; 1];
r22 = [-(L2-r2); 0; 0; 1];
gv = [0 -g 0 0];

G1 = -(m1 * gv * U11 * r11 + ...        % i=1, j=1
       m2 * gv * U21 * r22);            % i=1, j=2
G2 = -(m2 * gv * U22 * r22);            % i=2, j=2
G = simplify([G1; G2]);

%%
L2 = L1;

r1 = L1 / 2;
r2 = L2 / 2;
Iz1 = 1/3 * m1 * L1^2;
Iz2 = 1/3 * m2 * L2^2;

sM = simplify(eval(M));

sh1 = simplify(eval(h1));
sh2 = simplify(eval(h2));
sh = [sh1; sh2];

sg1 = simplify(eval(G1));
sg2 = simplify(eval(G2));
sg = [sg1; sg2];

%%
% [tau1; tau2] = M * [ddth1; ddth2] + h + G;
% [ddth1; ddth2] = inv(M) * ([tau1; tau2] - h - G);

syms tau1 tau2
DDTH = inv(M) * ([tau1; tau2] - h - G);
% DDTH1 = M\([tau1; tau2] - h - G);     -> same inverse matrix result

%%
dydt = simplify([dth1; DDTH(1); dth2; DDTH(2)]);
% matlabFunction(dydt, 'file', 'two_link.m', 'Optimize', false);

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
    tau1 = -dq1*1;
    tau2 = -dq2*1;
    
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




%% 1-DOF Simulation
% global m l g u;
% 
% dt = 0.005; ft = 5;
% q = pi/4; dq = 0;
% data = [];
% 
% m = 1; l = 1; g = 9.8148; n =1; u = 0;
% 
% FG = figure('Color', [1 1 1]);
% AX = axes('parent', FG);
% hold on;
% grid on;
% axis([-1.5 1.5 -1.5 1.5]);
% 
% Px = [0, 1];
% Py = [0, 0];
% 
% p = plot(Px, Py, '-ob', 'Linewidth', 3);
% 
% %% ODE45
% for cnt=0:dt:ft
%     [t,y] = ode45('one_link', [0 dt], [q; dq]);
%     
%     index = length(y);
%     
%     q = y(index, 1);
%     dq = y(index, 2);
%     
%     x = l*sin(q);
%     y = -l*cos(q);
%     Px = [0 x];
%     Py = [0 y];
%     
%     data(n,1) = cnt;
%     data(n,2) = q;
%     data(n,3) = dq;
%     n = n + 1;
%     
%     cmd = sprintf("Time: %2.2f", cnt);
%     clc
%     disp(cmd);
%     
%     if rem(n, 10) == 0
%         set(p, 'XData', Px, 'YData', Py);
%         drawnow
%     end
% end
% 
% %%
% FG2 = figure('Color', [1 1 1]);
% AX2 = axes('parent', FG2);
% grid on;
% 
% plot(data(:,1), data(:,2), 'r');
% hold on;
% plot(data(:,1), data(:,3), 'b');


















