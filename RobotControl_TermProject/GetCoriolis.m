%% Define
function H = GetCoriolis(q1, q2, dq1, dq2)
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

%% h term
% 2-DOF -> n=2

% Coriolis & Centrifugal term of first link
h111 = trace(U111 * J1 * U11.')...      % i=1, k=1, m=1, j=1
     + trace(U211 * J2 * U21.');        % i=1, k=1, m=1, j=2    dth1*dth1
h112 = trace(U212 * J2 * U21.');        % i=1, k=1, m=2, j=2    dth1*dth2
h121 = trace(U221 * J2 * U21.');        % i=1, k=2, m=1, j=2    dth2*dth1
h122 = trace(U222 * J2 * U21.');        % i=1, k=2, m=2, j=2    dth2*dth2

h1 = (dq1^2)*h111 + (dq1*dq2)*(h112 + h121) + (dq2^2)*h122;

% Coriolis & Centrifugal term of second link
h211 = trace(U211 * J2 * U22.');        % i=2, k=1, m=1, j=2
h212 = trace(U212 * J2 * U22.');        % i=2, k=1, m=2, j=2
h221 = trace(U221 * J2 * U22.');        % i=2, k=2, m=1, j=2
h222 = trace(U222 * J2 * U22.');        % i=2, k=2, m=2, j=2

h2 = (dq2^2)*h211 + (dq1*dq2)*(h212 + h221) + (dq2^2)*h222;

H = [h1; h2];
