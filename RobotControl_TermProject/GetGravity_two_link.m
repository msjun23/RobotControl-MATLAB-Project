%% Define
function G = GetGravity_two_link(q1, q2)
global L1 L2 m1 m2 r1 r2 g

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

%% Gravity term, C term
% syms r1 r2 Iz1 Iz2
r11 = [-(L1-r1); 0; 0; 1];
r22 = [-(L2-r2); 0; 0; 1];
gv = [0 -g 0 0];

G1 = -(m1 * gv * U11 * r11 + ...        % i=1, j=1
       m2 * gv * U21 * r22);            % i=1, j=2
G2 = -(m2 * gv * U22 * r22);            % i=2, j=2
G = [G1; G2];
