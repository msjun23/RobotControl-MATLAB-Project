%% Clear
clc
clear all
close all

%% DH Parameter, Tranformation matrix
syms L1 L2 L3 m1 m2 m3 r1 r2 r3 Iz1 Iz2 Iz3
syms th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3
syms g
syms tau1 tau2 tau3

dth = [dth1 dth2 dth3];         % Each joint angle

% Inertia of first body
I1xx = 0;
I1yy = Iz1;
I1zz = Iz1;

% Inertia of second body
I2xx = 0;
I2yy = Iz2;
I2zz = Iz2;

% Inertia of third body
I3xx = 0;
I3yy = Iz3;
I3zz = Iz3;

% 3-DOF DH parameter
d1 = 0;     d2 = 0;     d3 = 0;
a1 = L1;    a2 = L2;    a3 = L3;
al1 = 0;    al2 = 0;    al3 = 0;
% DH parameter matrix is as follows
% DH = [th1, d1, a1, al1;
%       th2, d2, a2, al2;
%       th3, d3, a3, al3];

% Function HT for get Homogeneous Transformation
T01 = HT(th1, d1, a1, al1);
T12 = HT(th2, d2, a2, al2);
T23 = HT(th3, d3, a3, al3);
T13 = T12 * T23;
T02 = T01 * T12;
T03 = T02 * T23;

%% Differential matrix(rotation)
Qr = [0 -1 0 0;
      1  0 0 0;
      0  0 0 0;
      0  0 0 0];

Q1 = Qr;
Q2 = Qr;
Q3 = Qr;

%% Velocity of a link, U matrix
U11 = Q1 * T01;             % i=1, j=1
U12 = zeros(4, 4);          % i=1, j=2
U13 = zeros(4, 4);          % i=1, j=3
U21 = Q1 * T02;             % i=2, j=1
U22 = T01 * Q2 * T12;       % i=2, j=2
U23 = zeros(4, 4);          % i=2, j=3
U31 = Q1 * T03;             % i=3, j=1
U32 = T01 * Q2 * T13;       % i=3, j=2
U33 = T02 * Q3 * T23;       % i=3, j=3

%% Kinetic energy of link i, Pseudo-inertia
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

J3(1, 1) = (-I3xx + I3yy + I3zz) / 2;
J3(2, 2) = ( I3xx - I3yy + I3zz) / 2;
J3(3, 3) = ( I3xx + I3yy - I3zz) / 2;
J3(1, 4) = -m3 * (L3 - r3);
J3(4, 1) = J3(1, 4);
J3(4, 4) = m3;

%% Inertia term, D term
% 3-DOF -> n=3
D11 = trace(U11 * J1 * U11.')...    % i=1, k=1, j=1
    + trace(U21 * J2 * U21.')...    % i=1, k=1, j=2
    + trace(U31 * J3 * U31.');      % i=1, k=1, j=3
M11 = simplify(D11);
D12 = trace(U22 * J2 * U21.')...    % i=1, k=2, j=2
    + trace(U32 * J3 * U31.');      % i=1, k=2, j=3
M12 = simplify(D12);
D13 = trace(U33 * J3 * U31.');      % i=1, k=3, j=3
M13 = simplify(D13);

M21 = M12;              % Inertia matrix -> symmetric
D22 = trace(U22 * J2 * U22.')...    % i=2, k=2, j=2
    + trace(U32 * J3 * U32.');      % i=2, k=2, j=3
M22 = simplify(D22);
D23 = trace(U33 * J3 * U32.');      % i=2, k=3, j=3
M23 = simplify(D23);

M31 = M13;              % Inertia matrix -> symmetric
M32 = M23;              % Inertia matrix -> symmetric
D33 = trace(U33 * J3 * U33.');      % i=3, k=3, j=3
M33 = simplify(D33);

M = [M11 M12 M13;
     M21 M22 M23;
     M31 M32 M33];

%% U matrices for get Coriolis & Centrifugal matrix
% 3-DOF -> n=3
U111 = Q1 * Q1 * T01;       % i=1, j=1, k=1
U112 = zeros(4, 4);         % i=1, j=1, k=2
U113 = zeros(4, 4);         % i=1, j=1, k=3
U121 = zeros(4, 4);         % i=1, j=2, k=1
U122 = zeros(4, 4);         % i=1, j=2, k=2
U123 = zeros(4, 4);         % i=1, j=2, k=3
U131 = zeros(4, 4);         % i=1, j=3, k=1
U132 = zeros(4, 4);         % i=1, j=3, k=2
U133 = zeros(4, 4);         % i=1, j=3, k=3

U211 = Q1 * Q1 * T02;       % i=2, j=1, k=1
U212 = Q1 * T01 * Q2 * T12; % i=2, j=1, k=2
U213 = zeros(4, 4);         % i=2, j=1, k=3
U221 = Q2 * T01 * Q2 * T12; % i=2, j=2, k=1
U222 = T01 * Q2 * Q2 * T12; % i=2, j=2, k=2
U223 = zeros(4, 4);         % i=2, j=2. k=3
U231 = zeros(4, 4);         % i=2, j=3, k=1
U232 = zeros(4, 4);         % i=2, j=3, k=2
U233 = zeros(4, 4);         % i=2, j=3, k=3

U311 = Q1 * Q1 * T03;       % i=3, j=1, k=1
U312 = Q1 * T01 * Q2 * T13; % i=3, j=1, k=2
U313 = Q1 * T02 * Q3 * T23; % i=3, j=1, k=3
U321 = Q1 * T01 * Q2 * T13;         % i=3, j=2, k=1
U322 = T01 * Q2 * Q2 * T13;         % i=3, j=2, k=2
U323 = T01 * Q2 * T12 * Q3 * T23;   % i=3, j=2, k=3
U331 = Q1 * T02 * Q3 * T23;         % i=3, j=3, k=1
U332 = T01 * Q2 * T12 * Q3 * T23;   % i=3, j=3, k=2
U333 = T02 * Q3 * Q3 * T23;         % i=3, j=3, k=3

%% Coriolis term, H term
% 3-DOF -> n=3
% Coriolis & Centrifugal term of first link
h111 = trace(U111 * J1 * U11.')...      % i=1, k=1, m=1, j=1
     + trace(U211 * J2 * U21.')...      % i=1, k=1, m=1, j=2
     + trace(U311 * J3 * U31.');        % i=1, k=1, m=1, j=3
h112 = trace(U212 * J2 * U21.')...      % i=1, k=1, m=2, j=2
     + trace(U312 * J3 * U31.');        % i=1, k=1, m=2, j=3
h113 = trace(U313 * J3 * U31.');        % i=1, k=1, m=3, j=3
h121 = trace(U221 * J2 * U21.')...      % i=1, k=2, m=1, j=2
     + trace(U321 * J3 * U31.');        % i=1, k=2, m=1, j=3
h122 = trace(U222 * J2 * U21.')...      % i=1, k=2, m=2, j=2
     + trace(U322 * J3 * U31.');        % i=1, k=2, m=2, j=3
h123 = trace(U323 * J3 * U31.');        % i=1, k=2, m=3, j=3
h131 = trace(U331 * J3 * U31.');        % i=1, k=3, m=1, j=3
h132 = trace(U332 * J3 * U31.');        % i=1, k=3, m=2, j=3
h133 = trace(U333 * J3 * U31.');        % i=1, k=3, m=3, j=3
h1_arr = [h111 h112 h113;
          h121 h122 h123;
          h131 h132 h133];
h1 = 0;
for k=1:3
    for m=1:3
        h1 = h1_arr(k,m) * dth(k) * dth(m) + h1;
    end
end

% Coriolis & Centrifugal term of second link
h211 = trace(U211 * J2 * U22.')...      % i=2, k=1, m=1, j=2
     + trace(U311 * J3 * U32.');        % i=2, k=1, m=1, j=3
h212 = trace(U212 * J2 * U22.')...      % i=2, k=1, m=2, j=2
     + trace(U312 * J3 * U32.');        % i=2, k=1, m=2, j=3
h213 = trace(U313 * J3 * U32.');        % i=2, k=1, m=3, j=3
h221 = trace(U221 * J2 * U22.')...      % i=2, k=2, m=1, j=2
     + trace(U321 * J3 * U32.');        % i=2, k=2, m=1, j=2
h222 = trace(U222 * J2 * U22.')...      % i=2, k=2, m=2, j=2
     + trace(U322 * J3 * U32.');        % i=2, k=2, m=2, j=3
h223 = trace(U323 * J3 * U32.');        % i=2, k=2, m=3, j=3
h231 = trace(U331 * J3 * U32.');        % i=2, k=3, m=1, j=3
h232 = trace(U332 * J3 * U32.');        % i=2, k=3, m=2, j=3
h233 = trace(U333 * J3 * U32.');        % i=2, k=3, m=3, j=3
h2_arr = [h211 h212 h213;
          h221 h222 h223;
          h231 h232 h233];
h2 = 0;
for k=1:3
    for m=1:3
        h2 = h2_arr(k,m) * dth(k) * dth(m) + h2;
    end
end

% Coriolis & Centrifugal term of third link
h311 = trace(U311 * J3 * U33.');        % i=3, k=1, m=1, j=3
h312 = trace(U312 * J3 * U33.');        % i=3, k=1, m=2, j=3
h313 = trace(U313 * J3 * U33.');        % i=3, k=1, m=3, j=3
h321 = trace(U321 * J2 * U33.');        % i=3, k=2, m=1, j=3
h322 = trace(U322 * J2 * U33.');        % i=3, k=2, m=2, j=3
h323 = trace(U323 * J3 * U33.');        % i=3, k=2, m=3, j=3
h331 = trace(U331 * J3 * U33.');        % i=3, k=3, m=1, j=3
h332 = trace(U332 * J3 * U33.');        % i=3, k=3, m=2, j=3
h333 = trace(U333 * J3 * U33.');        % i=3, k=3, m=3, j=3
h3_arr = [h311 h312 h313;
          h321 h322 h323;
          h331 h332 h333];
h3 = 0;
for k=1:3
    for m=1:3
        h3 = h3_arr(k,m) * dth(k) * dth(m) + h3;
    end
end

H = simplify([h1; h2; h3]);

%% Gravity term, G term
% Center of mass for each joint
r11 = [-(L1-r1); 0; 0; 1];
r22 = [-(L2-r2); 0; 0; 1];
r33 = [-(L3-r3); 0; 0; 1];

% Gravity vector
gv = [0 -g 0 0];

G1 = -(m1 * gv * U11 * r11 + ...        % i=1, j=1
       m2 * gv * U21 * r22 + ...        % i=1, j=2
       m3 * gv * U31 * r33);            % i=1, j=3
G2 = -(m2 * gv * U22 * r22 + ...        % i=2, j=2
       m3 * gv * U32 * r33);            % i=2, j=3
G3 = -(m3 * gv * U33 * r33);            % i=3, j=3

G = simplify([G1; G2; G3]);

%% Physical parameter
L3 = L2;
L2 = L1;
r1 = L1 / 2;
r2 = L2 / 2;
r3 = L3 / 2;
Iz1 = 1/3 * m1 * L1^2;
Iz2 = 1/3 * m2 * L2^2;
Iz3 = 1/3 * m3 * L3^2;

%% Dynamics Model
% [tau1; tau2; tau3] = M * [ddth1; ddth2; ddth3] + h + G;
% [ddth1; ddth2; tau3] = inv(M) * ([tau1; tau2; tau3] - h - G);

ddth = inv(M) * ([tau1; tau2; tau3] - H - G);

%% Define three_link.m function
dydt = simplify([dth1; ddth(1); dth2; ddth(2); dth3; ddth(3)]);
%matlabFunction(dydt, 'file', 'three_link.m', 'Optimize', false);
