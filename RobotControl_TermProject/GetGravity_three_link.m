%% Define
function G = GetGravity_three_link(q1, q2, q3)
global L1 L2 L3 m1 m2 m3 r1 r2 r3 g

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

%% Gravity term, C term
r11 = [-(L1-r1); 0; 0; 1];
r22 = [-(L2-r2); 0; 0; 1];
r33 = [-(L3-r3); 0; 0; 1];
gv = [0 -g 0 0];

G1 = -(m1 * gv * U11 * r11 + ...        % i=1, j=1
       m2 * gv * U21 * r22 + ...        % i=1, j=2
       m3 * gv * U31 * r33);            % i=1, j=3
G2 = -(m2 * gv * U22 * r22 + ...        % i=2, j=2
       m3 * gv * U32 * r33);            % i=2, j=3
G3 = -(m3 * gv * U33 * r33);            % i=3, j=3

G = [G1; G2; G3];
