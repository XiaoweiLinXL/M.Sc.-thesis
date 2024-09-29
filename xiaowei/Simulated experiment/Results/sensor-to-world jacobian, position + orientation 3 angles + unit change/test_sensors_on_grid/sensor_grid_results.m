%%
close all
clear
clc

%%
% Spheres
unit = 'decimeter';
if strcmp(unit, 'meter') 
    scale = 1;
elseif strcmp(unit, 'decimeter')
    scale = 0.1;
elseif strcmp(unit, 'centimeter')
    scale = 0.01;
elseif strcmp(unit, 'millimeter')
    scale = 0.001;
end

mu0 =  4*pi*1e-7; % air permeability
sph_dia = 3.175e-3*[0; 0; 1] / scale; % sphere diameter 3.175e-3 m (1/8 inch)

B_r = 13.2;
Volumn = (4/3)*pi*(norm(sph_dia)/2)^3; % volumn in the unit of (specified unit)^3

type = "3D";

%% Magnet configuration
magnet_conf = [0 0 0.2/scale 0 0 0;... 
               0 0 0.2/scale deg2rad(90) deg2rad(0) 0; ...
               0 0 0.2/scale deg2rad(0) deg2rad(90) 0].';

%% Sensor configuration
LB = [-0.2/scale, -0.2/scale];
UB = [0.2/scale, 0.2/scale];
z = 0;

number_of_sensor_on_side = 2:1:30;

sensor_grid_config = {};
for i=1:length(number_of_sensor_on_side)
    num_sample = number_of_sensor_on_side(i);
    
    C = cell(1, 2);
    [C{:}] = ndgrid(linspace(0, 1, num_sample));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    conf = [conf;zeros(1,size(conf,2));ones(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2))];

    sensor_grid_config{i} = conf;
end

% sensor_grid_config= [{[0 0.2/scale 0 1 0 0 0; -0.2/scale -0.2/scale 0 1 0 0 0; 0.2/scale -0.2/scale 0 1 0 0 0].'}, sensor_grid_config];
% sensor_grid_config = [{[-0.2/scale 0 0 1 0 0 0; 0.2/scale 0 0 1 0 0 0].'},sensor_grid_config];



%%
% Loop over each sensor configuration
reciprocal_condition_number_all_sensor_config = [];
min_singular_value_J_all_sensor_config = [];
min_singular_value_J_normalized_all_sensor_config = [];
min_singular_value_Jp_all_sensor_config = [];
min_singular_value_Jp_normalized_all_sensor_config = [];
min_singular_value_Jo_all_sensor_config = [];
min_singular_value_Jo_normalized_all_sensor_config = [];

for i=1:length(sensor_grid_config)
    
    
    % Extract one sensor configuration
    sensor_config = sensor_grid_config{i};
    sensor_number = size(sensor_config,2);
    sensor_position = sensor_config(1:3,:);
    sensor_orientation = sensor_config(4:7,:);
    magnitudes = vecnorm(sensor_orientation);
    sens_or_unitary = sensor_orientation ./ magnitudes;

    % Criteria
    reciprocal_condition_number = [];
    min_singular_value_J = [];
    min_singular_value_J_normalized = [];
    min_singular_value_Jp = [];
    min_singular_value_Jp_normalized = [];
    min_singular_value_Jo = [];
    min_singular_value_Jo_normalized = [];

    for j=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,j);
        theta = magnet_conf(4, j);
        phi = magnet_conf(5, j);
        psi = magnet_conf(6, j);

        J = [];
        for k = 1:sensor_number
            J = [J;J_analytical_sensor_B_to_magnet_conf_euler_angle_XYX_scaled_uni(sensor_position(:,k), sensor_orientation(:,k), ...
                magnet_pos, B_r, Volumn, theta, phi, psi, type)];
        end

        Jp = J(:,1:3);
        Jo = J(:,4:6);

        num_dof = 5;
        num_dof_position = 3;
        num_dof_rotation = 2;

        sigma_J = svd(J);
        sigma_Jp = svd(Jp);
        sigma_Jo = svd(Jo);

        rcond = sigma_J(num_dof)/sigma_J(1);
        min_svd_J = sigma_J(num_dof);
        min_svd_J_normalized = min_svd_J/sqrt(sensor_number);
        min_svd_Jp = sigma_Jp(num_dof_position);
        min_svd_Jp_normalized = min_svd_Jp/sqrt(sensor_number);
        min_svd_Jo = sigma_Jo(num_dof_rotation);
        min_svd_Jo_normalized  = min_svd_Jo/sqrt(sensor_number);

        reciprocal_condition_number = [reciprocal_condition_number; rcond];
        min_singular_value_J = [min_singular_value_J; min_svd_J];
        min_singular_value_J_normalized = [min_singular_value_J_normalized; min_svd_J_normalized];
        min_singular_value_Jp = [min_singular_value_Jp; min_svd_Jp];
        min_singular_value_Jp_normalized = [min_singular_value_Jp_normalized; min_svd_Jp_normalized];
        min_singular_value_Jo = [min_singular_value_Jo; min_svd_Jo];
        min_singular_value_Jo_normalized = [min_singular_value_Jo_normalized; min_svd_Jo_normalized];
    end

    % Put the minimum among all magnet configs into the sensor config
    % collection
    reciprocal_condition_number_all_sensor_config = [reciprocal_condition_number_all_sensor_config; min(reciprocal_condition_number)];
    min_singular_value_J_all_sensor_config = [min_singular_value_J_all_sensor_config; min(min_singular_value_J)];
    min_singular_value_J_normalized_all_sensor_config = [min_singular_value_J_normalized_all_sensor_config; min(min_singular_value_J_normalized)];
    min_singular_value_Jp_all_sensor_config = [min_singular_value_Jp_all_sensor_config; min(min_singular_value_Jp)];
    min_singular_value_Jp_normalized_all_sensor_config = [min_singular_value_Jp_normalized_all_sensor_config; min(min_singular_value_Jp_normalized)];
    min_singular_value_Jo_all_sensor_config = [min_singular_value_Jo_all_sensor_config; min(min_singular_value_Jo)];
    min_singular_value_Jo_normalized_all_sensor_config = [min_singular_value_Jo_normalized_all_sensor_config; min(min_singular_value_Jo_normalized)];
end

%%
number_of_sensor = [number_of_sensor_on_side.^2];

figure
plot(number_of_sensor, reciprocal_condition_number_all_sensor_config)
xlabel('number of sensors')
ylabel('min reciprocal condition number')
ylim([0,1])
title('Min reciprocal condition number vs number of sensor on one row')
grid on


figure
plot(number_of_sensor, min_singular_value_J_all_sensor_config)
xlabel('number of sensors')
ylabel('min min singular value of J')
title('Min min singular value of J vs number of sensor on one row')
grid on

figure
plot(number_of_sensor, min_singular_value_J_normalized_all_sensor_config)
xlabel('number of sensors')
ylabel('min min singular value of J')
title('Min min singular value of J normalized vs number of sensor on one row')
grid on


figure
plot(number_of_sensor, min_singular_value_Jp_all_sensor_config)
xlabel('number of sensors')
ylabel('min min singular value of Jp')
title('Min min singular value of Jp vs number of sensor on one row')
grid on

figure
plot(number_of_sensor, min_singular_value_Jp_normalized_all_sensor_config)
xlabel('number of sensors')
ylabel('min min singular value of Jp')
title('Min min singular value of Jp normalized vs number of sensor on one row')
grid on


figure
plot(number_of_sensor, min_singular_value_Jo_all_sensor_config)
xlabel('number of sensors')
ylabel('min min singular value of Jo')
title('Min min singular value of Jo vs number of sensor on one row')
grid on

figure
plot(number_of_sensor, min_singular_value_Jo_normalized_all_sensor_config)
xlabel('number of sensors on one row')
ylabel('min min singular value of Jo')
title('Min min singular value of Jo normalized vs number of sensor on one row')
grid on




%% Functions
% Quaternion to rotation matrix
function R = quaternionToMatrix(q)
    % Extract the scalar and vector parts from the quaternion
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);

    % Calculate the elements of the rotation matrix
    R11 = 1 - 2*y^2 - 2*z^2;
    R12 = 2*x*y - 2*z*w;
    R13 = 2*x*z + 2*y*w;
    R21 = 2*x*y + 2*z*w;
    R22 = 1 - 2*x^2 - 2*z^2;
    R23 = 2*y*z - 2*x*w;
    R31 = 2*x*z - 2*y*w;
    R32 = 2*y*z + 2*x*w;
    R33 = 1 - 2*x^2 - 2*y^2;

    % Combine the elements into the rotation matrix
    R = [R11, R12, R13;
         R21, R22, R23;
         R31, R32, R33];
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame, euler angle
function J = J_analytical_sensor_B_to_magnet_conf_euler_angle_XYX_scaled_unit(sens_pos, sens_or, magnet_pos, B_r, Volumn, theta, phi, psi, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    Rx = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];
    
    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];

    Rx_2 = [1 0 0;                  % rotate about x-axis
    0 cos(psi) -sin(psi);
    0 sin(psi) cos(psi)];

    mu_world_hat = Rx * Ry * Rx_2 * [0;0;1];

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r_world = sens_pos - magnet_pos;
    r_world_hat = r_world/norm(r_world);
    
    J_position = -((3*B_r*Volumn)/(4*pi*norm(r_world)^4)) * ((eye(3)-5*(r_world_hat*r_world_hat.'))*(r_world_hat.'*mu_world_hat) ...
                + mu_world_hat*r_world_hat.' + r_world_hat*mu_world_hat.');

    % Second part of the jacobian, from sensor reading to magnet
    % orientation

    % Extract rotation information
    dmu_dangle = [0 cos(phi)*cos(psi) -sin(phi)*sin(psi); 
                          sin(theta)*sin(psi)-cos(phi)*cos(psi)*cos(theta) sin(phi)*cos(psi)*sin(theta) -cos(theta)*cos(psi)+cos(phi)*sin(psi)*sin(theta); 
                          -cos(phi)*cos(psi)*sin(theta)-sin(psi)*cos(theta) -sin(phi)*cos(psi)*cos(theta) -cos(phi)*sin(psi)*cos(theta)-cos(psi)*sin(theta)];

    J_angle = (B_r*Volumn/(4*pi))*(1/(norm(r_world)^3))*(3*(r_world_hat*r_world_hat.')-eye(3))*dmu_dangle;

    J = [J_position, J_angle];

    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end



