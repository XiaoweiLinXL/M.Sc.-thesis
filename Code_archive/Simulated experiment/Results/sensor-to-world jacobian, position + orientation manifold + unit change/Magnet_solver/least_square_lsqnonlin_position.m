close all
clear
clc

%% lsqnonlin
sph_dia = 3.175e-3*[0; 0; 1];
Volumn = (4/3)*pi*(norm(sph_dia)/2)^3;
B_r = 1.32;
p_star = [0;0;0.05];
R_star = eye(3);
sens_pos_collection = [-0.28,-0.28,0,-0.28,0.28,0,0.28,-0.28,0,0.28,0.28,0,4];
B_star = mag_field_vector(sens_pos_collection, p_star, B_r, Volumn, R_star);

x0 = [0;0;0.2];

options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt', ...
'FunctionTolerance', 1e-30, 'OptimalityTolerance', 1e-30, 'StepTolerance', ...
1e-30, 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6);
fun = @(x) residual_fun(x, sens_pos_collection, B_star, B_r, Volumn);
[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[],[],options);

%% Residual function
function residual = residual_fun(magnet_pos, sens_pos_collection, B_star, B_r, Volumn)
    B_estimation = mag_field_vector(sens_pos_collection, magnet_pos, B_r, Volumn, eye(3));
    residual = B_star - B_estimation;
end


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