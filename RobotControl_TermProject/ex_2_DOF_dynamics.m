%% Clear
clc
clear all
close all

%% DH Parameter, Tranformation matrix
syms L1 L2 m1 m2 r1 r2 Iz1 Iz2
syms th1 th2 dth1 dth2 ddth1 ddth2

% Inertia of first body
I1xx = 0;
I1yy = Iz1;
I1zz = Iz1;

% Inertia of second body
I2xx = 0;
I2yy = Iz2;
I2zz = Iz2;

% 2-DOF DH parameter
d1 = 0;     d2 = 0;
a1 = L1;    a2 = L2;
al1 = 0;    al2 = 0;
DH = [th1, d1, a1, al1 ;
      th2, d2, a2, al2];

% Function HT for get Homogeneous Transformation
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

% dT01t1 = Qr * T01;
% dT01t2 = diff(T01);       <- same result

%% Velocity of a link, U matrix
U11 = Q1 * T01;             % i=1, j=1
U12 = zeros(4, 4);          % i=1, j=2
U21 = Q1 * T02;             % i=2, j=1
U22 = T01 * Q2 * T12;       % i=2, j=2

%% Kinetic energy of link i, Pseudo-inverse
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

%% Inertia matrix, D term
% 2-DOF -> n=2
D11 = trace(U11 * J1 * U11.')...    % i=1, k=1, j=1
    + trace(U21 * J2 * U21.');      % i=1, k=1, j=2
M11 = simplify(D11);

D12 = trace(U22 * J2 * U21.');      % i=1, k=2, j=2
M12 = simplify(D12);

% Inertia matrix -> symmetric
M21 = M12;

D22 = trace(U22 * J2 * U22.');      % i=2, k=2, j=2
M22 = simplify(D22);

M = [M11 M12; M21 M22];

%% Inertia matrix -> Define using summation function
% n = 2;
% for i = 1:n
%     for k = 1:n
%         nM(i,k) = Inertia(i, k, n, U11, U12, U21, U22, J1, J2);
%     end
% end

%% Coriolis & Centrifugal matrix
% 2-DOF -> n=2
U111 = Q1 * Q1 * T01;       % i=1, j=1, k=1
U112 = zeros(4, 4);         % i=1, j=1, k=2
U121 = zeros(4, 4);         % i=1, j=2, k=1
U122 = zeros(4, 4);         % i=1, j=2, k=2
U211 = Q1 * Q1 * T02;       % i=2, j=1, k=1
U212 = Q1 * T01 * Q2 * T12; % i=2, j=1, k=2
U221 = Q2 * T01 * Q2 * T12; % i=2, j=2, k=1
U222 = T01 * Q2 * Q2 * T12; % i=2, j=2, k=2

%%  Coriolis & Centrifugal matrix -> Define using function
% n = 2;
% 
% for i = 1:n
%     for j = 1:n
%         for k = 1:n
%             cmd = sprintf("nU%d%d%d = dUdq(i, j, k, T01, T12, T02, Q1, Q2);", i, j, k);
%             eval(cmd);
%         end
%     end
% end

%% h term
% 2-DOF -> n=2

% Coriolis & Centrifugal term of first link
h111 = trace(U111 * J1 * U11.')...      % i=1, k=1, m=1, j=1
     + trace(U211 * J2 * U21.');        % i=1, k=1, m=1, j=2    dth1*dth1
h112 = trace(U212 * J2 * U21.');        % i=1, k=1, m=2, j=2    dth1*dth2
h121 = trace(U221 * J2 * U21.');        % i=1, k=2, m=1, j=2    dth2*dth1
h122 = trace(U222 * J2 * U21.');        % i=1, k=2, m=2, j=2    dth2*dth2

h1 = (dth1^2)*h111 + (dth1*dth2)*(h112 + h121) + (dth2^2)*h122;

% Coriolis & Centrifugal term of second link
h211 = trace(U211 * J2 * U22.');        % i=2, k=1, m=1, j=2
h212 = trace(U212 * J2 * U22.');        % i=2, k=1, m=2, j=2
h221 = trace(U221 * J2 * U22.');        % i=2, k=2, m=1, j=2
h222 = trace(U222 * J2 * U22.');        % i=2, k=2, m=2, j=2

h2 = (dth2^2)*h211 + (dth1*dth2)*(h212 + h221) + (dth2^2)*h222;

h = simplify([h1; h2]);

%% Gravity term, C term
% syms r1 r2 Iz1 Iz2
syms g
r11 = [-(L1-r1); 0; 0; 1];
r22 = [-(L2-r2); 0; 0; 1];
gv = [0 -g 0 0];

G1 = -(m1 * gv * U11 * r11 + ...        % i=1, j=1
       m2 * gv * U21 * r22);            % i=1, j=2
G2 = -(m2 * gv * U22 * r22);            % i=2, j=2
G = simplify([G1; G2]);

%% Physical parameter
L2 = L1;
r1 = L1 / 2;
r2 = L2 / 2;
Iz1 = 1/3 * m1 * L1^2;
Iz2 = 1/3 * m2 * L2^2;

% Inertia term
sM = simplify(eval(M));

% Coriolis & Centrifugal term
sh1 = simplify(eval(h1));
sh2 = simplify(eval(h2));
sh = [sh1; sh2];

% Gravity term
sg1 = simplify(eval(G1));
sg2 = simplify(eval(G2));
sg = [sg1; sg2];

%%
% [tau1; tau2] = M * [ddth1; ddth2] + h + G;
% [ddth1; ddth2] = inv(M) * ([tau1; tau2] - h - G);

syms tau1 tau2
DDTH = inv(M) * ([tau1; tau2] - h - G);
% DDTH = M\([tau1; tau2] - h - G);     -> same inverse matrix result

%% Define two_link.m function
dydt = simplify([dth1; DDTH(1); dth2; DDTH(2)]);
%matlabFunction(dydt, 'file', 'two_link.m', 'Optimize', false);
