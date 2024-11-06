%%
close all
clear
clc

%%
syms theta phi psi

Rz = [cos(theta) -sin(theta) 0; 
      sin(theta) cos(theta) 0;
      0 0 1];

Rx_1 = [1 0 0;
        0 cos(theta) -sin(theta);
        0 sin(theta) cos(theta)];

Ry = [cos(phi) 0 sin(phi);
      0 1 0;
      -sin(phi) 0 cos(phi)];

Rx = [1 0 0;
      0 cos(psi) -sin(psi);
      0 sin(psi) cos(psi)];

mu_hat_zyx = Rz*Ry*Rx*[0;0;1];

Rotation_xyx = Rx_1*Ry*Rx;

mu_hat_xyx = Rotation_xyx*[0;0;1];

J_zyx = jacobian(mu_hat_zyx, [theta, phi, psi]);

J_zyx_sub = subs(J_zyx, [theta, phi, psi], [0,pi/2,pi/2]);

mu_hat_zyx_sub = subs(mu_hat_zyx, [theta, phi, psi], [0,pi/2,pi/2]);

J_xyx = jacobian(mu_hat_xyx, [theta, phi, psi]);

J_xyx_sub = subs(J_xyx, [theta, phi, psi], [0,0,0]);