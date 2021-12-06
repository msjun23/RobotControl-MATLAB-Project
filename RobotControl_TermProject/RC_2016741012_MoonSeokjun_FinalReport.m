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















