close all
clear 
clc

%%
syms w x y z real

q = [w,x,y,z];

Rotation = quaternionToMatrix(q);

v = [0;0;1];

Cv = Rotation*[0;0;1];

J = jacobian(Cv, [w,x,y,z]);
%%
quaternionToMatrix([0,0,1,0])


%% General function
% Quaternion to rotation matrix
function R = quaternionToMatrix(q)
    % Extract the scalar and vector parts from the quaternion
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);

    % Calculate the elements of the rotation matrix
    R11 = w^2 + x^2 - y^2 - z^2;
    R12 = 2*x*y - 2*z*w;
    R13 = 2*x*z + 2*y*w;
    R21 = 2*x*y + 2*z*w;
    R22 = w^2 + y^2 - x^2 - z^2;
    R23 = 2*y*z - 2*x*w;
    R31 = 2*x*z - 2*y*w;
    R32 = 2*y*z + 2*x*w;
    R33 = w^2 + z^2 - x^2 - y^2;

    % Combine the elements into the rotation matrix
    R = [R11, R12, R13;
         R21, R22, R23;
         R31, R32, R33];
end