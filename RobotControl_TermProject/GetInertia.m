%% Define
function D = GetInertia(q1, q2)
global L1 L2 m1 m2 r1 r2 Iz1 Iz2

%% DH Parameter, Tranformation matrix
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
% DH = [q1, d1, a1, al1 ;
%       q2, d2, a2, al2];

% Function HT for get Homogeneous Transformation
T01 = HT(q1, d1, a1, al1);
T12 = HT(q2, d2, a2, al2);
T02 = T01 * T12;

%% Differential matrix(rotation)
Qr = [0 -1 0 0;
      1  0 0 0;
      0  0 0 0;
      0  0 0 0];

Q1 = Qr;
Q2 = Qr;

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

D12 = trace(U22 * J2 * U21.');      % i=1, k=2, j=2

% Inertia matrix -> symmetric
D21 = D12;

D22 = trace(U22 * J2 * U22.');      % i=2, k=2, j=2

D = [D11 D12; D21 D22];
