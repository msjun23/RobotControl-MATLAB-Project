%% Define
function D = GetInertia_three_link(q1, q2, q3)
global L1 L2 L3 m1 m2 m3 r1 r2 r3 Iz1 Iz2 Iz3

%% DH Parameter, Tranformation matrix
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

% Function HT for get Homogeneous Transformation
T01 = HT(q1, d1, a1, al1);
T12 = HT(q2, d2, a2, al2);
T23 = HT(q3, d3, a3, al3);
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

J3(1, 1) = (-I3xx + I3yy + I3zz) / 2;
J3(2, 2) = ( I3xx - I3yy + I3zz) / 2;
J3(3, 3) = ( I3xx + I3yy - I3zz) / 2;
J3(1, 4) = -m3 * (L3 - r3);
J3(4, 1) = J3(1, 4);
J3(4, 4) = m3;

%% Inertia matrix, D term
% 3-DOF -> n=3
D11 = trace(U11 * J1 * U11.')...    % i=1, k=1, j=1
    + trace(U21 * J2 * U21.')...    % i=1, k=1, j=2
    + trace(U31 * J3 * U31.');      % i=1, k=1, j=3
D12 = trace(U22 * J2 * U21.')...    % i=1, k=2, j=2
    + trace(U32 * J3 * U31.');      % i=1, k=2, j=3
D13 = trace(U33 * J3 * U31.');      % i=1, k=3, j=3

D21 = D12;              % Inertia matrix -> symmetric
D22 = trace(U22 * J2 * U22.')...    % i=2, k=2, j=2
    + trace(U32 * J3 * U32.');      % i=2, k=2, j=3
D23 = trace(U33 * J3 * U32.');      % i=2, k=3, j=3

D31 = D13;              % Inertia matrix -> symmetric
D32 = D23;              % Inertia matrix -> symmetric
D33 = trace(U33 * J3 * U33.');      % i=3, k=3, j=3

D = [D11 D12 D13;
     D21 D22 D23;
     D31 D32 D33];
 