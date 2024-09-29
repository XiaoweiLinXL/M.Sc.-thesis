close all 
clear
clc

%% 
syms w1 w2 w3 real

w = [w1;w2;w3];

w_norm = sqrt(w1^2+w2^2+w3^2);

normalizd_w = w/w_norm;

Rotation = eye(3) + (sin(w_norm)/w_norm)*skew_symmetric(w) + ((1-cos(w_norm))/w_norm^2)*skew_symmetric(w)*skew_symmetric(w);

Rotation_Rodrigues = eye(3) + sin(w_norm)*skew_symmetric(normalizd_w) + ...
(1-cos(w_norm))*skew_symmetric(normalizd_w)*skew_symmetric(normalizd_w);

v = [0;0;1];

Cv_Rodrigues = Rotation_Rodrigues*v;

J_Rodrigues = jacobian(Cv_Rodrigues, [w1,w2,w3]);

J_Rodrigues_sub = subs(J_Rodrigues, [w1,w2,w3],[2*pi,0,0]);

w_bar = [pi/3,pi/6,pi/2];

R0 = subs(Rotation,[w1,w2,w3],w_bar);

G1 = skew_symmetric([1;0;0]);
G2 = skew_symmetric([0;1;0]);
G3 = skew_symmetric([0;0;1]);

J=[G1*R0*v G2*R0*v G3*R0*v];

dRp_dR = -R0*skew_symmetric(v);

% Cv = Rotation*v;
% 
% J = jacobian(Cv, [w1,w2,w3]);
% 
% J_sub = subs(J,[w1,w2,w3],w_bar);

v_cross = skew_symmetric(R0*v);

J1 = -Rotation*v_cross;

J1_sub = subs(J1, [w1,w2,w3],w_bar);

% J2 = eye(3)-((1-cos(w_norm))/(w_norm^2))*skew_symmetric(w)+((w_norm-sin(w_norm))/(w_norm^3))*skew_symmetric(w)^2;







%% Function
function x_cross = skew_symmetric(x)
    x_cross = [0 -x(3) x(2);
               x(3) 0 -x(1);
               -x(2) x(1) 0];
end