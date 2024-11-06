%%
close all
clear
clc

%%
syms theta phi psi

Rx_1 = [1 0 0;
        0 cos(theta) -sin(theta); 
        0 sin(theta) cos(theta)];

Ry = [cos(phi) 0 sin(phi);
      0 1 0;
      -sin(phi) 0 cos(phi)];

Rx_2 = [1 0 0;
        0 cos(psi) -sin(psi);
        0 sin(psi) cos(psi)];

Rotation = Rx_2 * Ry * Rx_1;

v = [0;0;1];

Cv = Rotation * v;

J_xyx_extrinsic = jacobian(Cv, [theta, phi, psi]);

skew_Cv_extrinsic = [0 -Cv(3) Cv(2); Cv(3) 0 -Cv(1); -Cv(2) Cv(1) 0];

J_xyx_extrinsic_analytic = skew_Cv_extrinsic*[Rx_2*Ry*[1;0;0] Rx_2*[0;1;0] [1;0;0]];
simplify(J_xyx_extrinsic_analytic)
