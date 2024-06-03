syms theta phi psi


Rx_1 = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Ry = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
Rx_2 = [1 0 0; 0 cos(psi) -sin(psi); 0 sin(psi) cos(psi)];

Rotation = Rx_1 * Ry * Rx_2;

mu_norm = 0.176;
mu0 = 4*pi*1e-7;
r1 = [0.6;0.2;-0.65];
r1_hat = r1/norm(r1);

B1 = (3*r1_hat*r1_hat.' - eye(3)) * Rotation * [0;0;1];

J1 = jacobian(B1, [theta, phi, psi]);

J1 = subs(J1, [theta,phi,psi],[0,pi/3,0])

r2 = [0;0;0.1];
r2_hat = r2/norm(r2);

B2 = (3*r2_hat*r2_hat.' - eye(3)) * Rotation * [0;0;1];

J2 = jacobian(B2, [theta, phi, psi])

J2 = subs(J2, [theta,phi,psi],[0,pi/3,0])

r3 = [-0.4;-3.7;-1.5];
r3_hat = r3/norm(r3);

B3 = (3*r3_hat*r3_hat.' - eye(3)) * Rotation * [0;0;1];

J3 = jacobian(B3, [theta, phi, psi])

J3 = subs(J3, [theta,phi,psi],[0,pi/3,0])

r4 = [4.9;9.6;-5.5];
r4_hat = r4/norm(r4);

B4 = (3*r4_hat*r4_hat.' - eye(3)) * Rotation * [0;0;1];

J4 = jacobian(B4, [theta, phi, psi])

J4 = subs(J4, [theta,phi,psi],[0,pi/3,0])

r5 = [-7.6;12.9;2.7];
r5_hat = r5/norm(r5);

B5 = (3*r5_hat*r5_hat.' - eye(3)) * Rotation * [0;0;1];

J5 = jacobian(B5, [theta, phi, psi])

J5 = subs(J5, [theta,phi,psi],[0,pi/3,0])

% subs(Rotation, [theta,phi,psi],[0,pi/3,0])
% 
% subs(Rotation, [theta,phi,psi],[pi/4,0,pi/4])
