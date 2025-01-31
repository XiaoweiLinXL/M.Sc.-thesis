close all
clear 
clc

%% Constant
scale = 1;
sph_dia = 3.175e-3*[0; 0; 1]/scale;
Volumn = 24*(4/3)*pi*(norm(sph_dia)/2)^3;
B_r = 1.32;

type = "3D";

plot_progress = false;
iteration = 1000;

% sens_pos_collection = [0.025, 0,0,-0.025,0,0,0,0.025,0,0,-0.025,0]/scale;
sens_pos_collection = [0.05, 0,0,-0.05,0,0,0,0.05,0,0,-0.05,0]/scale;
% sens_pos_collection = [0.1120, 0,0,-0.1120,0,0,0,0.1120,0,0,-0.1120,0]/scale;
% sens_pos_collection = [0.1250, 0,0,-0.1250,0,0,0,0.1250,0,0,-0.1250,0]/scale;
% sens_pos_collection = [0.2150, 0,0,-0.2150,0,0,0,0.2150,0,0,-0.2150,0]/scale;
% sens_pos_collection = [0.3550, 0,0,-0.3550,0,0,0,0.3550,0,0,-0.3550,0]/scale;

sens_pos_collection = [sens_pos_collection, 4];
magnet_pos = [0;0;0.15];
R_star = rotx(45);


B1 = mag_field_vector(sens_pos_collection, magnet_pos, B_r, Volumn, R_star)

magnet_pos = magnet_pos + [0.0;0.0;0.0];
R_star =  R_star * rotx(45);
B2 = mag_field_vector(sens_pos_collection, magnet_pos, B_r, Volumn, R_star)

norm(B1-B2)

%% Magnetic field vector for all sensors
function B_vector = mag_field_vector(sens_pos_collection, magnet_pos, B_r, Volumn, R_star)
    B_vector = [];
    sens_num = sens_pos_collection(end);
    sens_pos_collection(end) = [];
    sens_pos_collection  = reshape(sens_pos_collection, 3, []);
    for sens_index = 1:sens_num
        sens_pos = sens_pos_collection(:,sens_index);
        B_vector = [B_vector; mag_field(sens_pos, magnet_pos, B_r, Volumn, R_star)];
    end
end

%% Magnetic field
function B = mag_field(sens_pos, magnet_pos, B_r, Volumn, R_star)
    r = sens_pos-magnet_pos;
    r_hat = r/norm(r);
    B = (B_r*Volumn)/(4*pi*(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*R_star*[0;0;1];
end