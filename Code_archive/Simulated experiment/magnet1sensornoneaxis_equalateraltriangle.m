%%
close all
clear
clc
%%
% Spheres
mu0 =  4*pi*1e-7; % air permeability
sph_dia = 3.175e-3*[0; 0; 1]; % sphere diameter
sph_dip = 13.2*4/3*pi*(norm(sph_dia)/2)^3/mu0; % spheres magnetic dipole
mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole

delta = 1e-7;
type = "1D";

%% Magnet configurations
% Workspace as a box in m
LB = [-0.05, -0.05, 0.3];
UB = [0.05, 0.05, 0.4];

% Number of samples within the box will be num_samp^3
num_samp = 0;

C = cell(1, 3);
[C{:}] = ndgrid(linspace(0, 1, num_samp));
C = cellfun(@(a) a(:), C, 'Uni', 0);
combo = [C{:}];
conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
    repmat(LB, [size(combo, 1), 1]);
conf = conf.';

magnet_conf = conf;

magnet_conf = [[0;0;0.2059]];
%% Genetic algorithm
for i = 1:1:1
    lb = [0.2 0.0 0.0 -1 -1 -1 -1 ...
          -0.2*0.5 0.2*sqrt(3)/2 0.0 -1 -1 -1 -1 ...
          -0.2*0.5 -0.2*sqrt(3)/2 0.0 -1 -1 -1 -1 ...
          -0.2 -0.2 0.0 -1 -1 -1 -1 ...
          -0.2 -0.2 0.0 -1 -1 -1 -1 ...
          -0.2 -0.2 0.0 -1 -1 -1 -1 ...
          -0.2 -0.2 0.0 -1 -1 -1 -1 ...
          -0.2 -0.2 0.0 -1 -1 -1 -1 ...
          -0.2 -0.2 0.0 -1 -1 -1 -1 ...
          -0.2 -0.2 0.0 -1 -1 -1 -1 ...
          3];
    ub = [0.2 0.0 0.0 1 1 1 1 ...
          -0.2*0.5 0.2*sqrt(3)/2 0.0 1 1 1 1 ...
          -0.2*0.5 -0.2*sqrt(3)/2 0.0 1 1 1 1 ...
          0.2 0.2 0.0 1 1 1 1 ...
          0.2 0.2 0.0 1 1 1 1 ...
          0.2 0.2 0.0 1 1 1 1 ...
          0.2 0.2 0.0 1 1 1 1 ...
          0.2 0.2 0.0 1 1 1 1 ...
          0.2 0.2 0.0 1 1 1 1 ...
          0.2 0.2 0.0 1 1 1 1 ...
          3];
    options = optimoptions(@gamultiobj,'Display','iter','PopulationSize',500,'MaxStallGeneration',500);
    
    fun = @(x) min_fun(x, magnet_conf, mu, type);
    [sol, fval, exitflag, output] = gamultiobj(fun, length(lb), [], [], [], [], lb, ub, [], length(lb), options);
    
    filename = sprintf('data_%d.mat', i); % Generate file name
    save(filename); % Save data, replace 'data_to_save' with your data variable

    
    % Evaluate
    load(filename)
    sens_conf = [sol];
    [obj_rcond_all_solution, obj_B_all_solution, min_rcond, min_B, min_singular_value] = evaluate(sens_conf, magnet_conf, mu, type);
    sol
    obj_rcond_all_solution
    obj_B_all_solution
    min_rcond
    min_B
    min_singular_value
end

%%
load("data_1.mat")
sens_conf = [sol];
[obj_rcond_all_solution, obj_B_all_solution, min_rcond, min_B, min_singular_value] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond_all_solution
obj_B_all_solution
min_rcond
min_B
min_singular_value

%% Plot the pareto plane 2D
optimal_min_rcond = min(obj_rcond_all_solution, [], 1);
optimal_min_B = min(obj_B_all_solution, [], 1);
pareto_plane = [optimal_min_rcond;optimal_min_B];
figure
plot(pareto_plane(1,:), pareto_plane(2,:), 'bo', 'MarkerSize', 10);
xlabel("Min rcond in workspace")
ylabel("Min B in workspace(Tesla)")
title("Pareto plane from the solver")
grid on

%% Plot the pareto plane 3D
optimal_sensor_num = sol(:,end);
optimal_min_rcond = min(obj_rcond_all_solution, [], 1);
optimal_min_B = min(obj_B_all_solution, [], 1);
pareto_plane = [optimal_min_rcond;optimal_min_B;optimal_sensor_num'];
figure
plot3(pareto_plane(1,:), pareto_plane(2,:), pareto_plane(3,:), 'b.', 'MarkerSize', 10);
xlabel("Min rcond in workspace")
ylabel("Min B in workspace(Tesla)")
zlabel("# of sensor")
title("Pareto plane from the solver")
grid on

%% Jacobian test
theta = pi/6.3;
phi = 0;
gamma = 0;
Rx = [1 0 0;                  % rotate about x-axis
0 cos(theta) -sin(theta);
0 sin(theta) cos(theta)];

Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
0 1 0;
-sin(phi) 0 cos(phi)];

Rz = [cos(gamma) -sin(gamma) 0;
      sin(gamma) cos(gamma) 0;
      0 0 1];

sensor_output([0.0;0.0;0.0], [1;1;1;1], [0.1;0.3;0.3], Rx*Ry*Rz*[0;0;0.0176], '3D')
sensor_output_from_angle([0.0;0.0;0.0], [1;1;1;1], [0.1;0.3;0.3], 0.0176, theta, phi, '3D')
sensor_output_from_quaternion([0.0;0.0;0.0], [1;1;1;1], [0.1;0.3;0.3], 0.0176, cos(theta/2), sin(theta/2), 0, 0, '3D')

jacobian_numerical_sensor_reading_to_magnet_conf_euler_angle([-0.03;0.45;0.23], [1;0;0;0], [0;0;0.3], 0.0176, theta, phi, 1e-7, '3D')
jacobian_B_sensor_to_magnet_conf([-0.03;0.45;0.23], [1;0;0;0], [0;0;0.3], Rx*Ry*Rz*[0;0;0.0176], theta, phi, '3D')

J_numerical_1 = jacobian_numerical_sensor_reading_to_magnet_conf_quaternion([-0.03;0.45;0.23], [0.5, 0.5, 0.5, 0.5,], [0;0;0.3], 0.0176, 0.5, 0.5, 0.5, 0.5, 1e-7, '3D')
J_analytical_1 = jacobian_B_sensor_to_magnet_conf_quaternion([3;0.1;7], [0.5, 0.5, 0.5, 0.5,], [0;0;0.3], quaternionToMatrix([0.5, 0.5, 0.5, 0.5,])*[0;0;0.0176], 0.5, 0.5, 0.5, 0.5, '3D')
J_numerical_2 = jacobian_numerical_sensor_reading_to_magnet_conf_quaternion([0.9;6.1;3.1], [0.5, 0.5, 0.5, 0.5,], [0;0;0.3], 0.0176, 0.5, 0.5, 0.5, 0.5, 1e-7, '3D')
J_analytical_2 = jacobian_B_sensor_to_magnet_conf_quaternion([0.9;6.1;3.1], [0.5, 0.5, 0.5, 0.5,], [0;0;0.3], quaternionToMatrix([0.5, 0.5, 0.5, 0.5,])*[0;0;0.0176], 0.5, 0.5, 0.5, 0.5, '3D')
%% General function
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

% Magnetic field
function B = mag_field(r, mu)
% Compute magnetic field at location p for a magnet of magnetic dipole
% mu
    mu0 = 4*pi*1e-7;
    hatp = r/norm(r);

    B = mu0/(4*pi*norm(r)^3)*(3*(hatp*hatp.') - eye(3))*mu;
end

% Sensor output based on sensor type
function output = sensor_output(sens_pos, sens_or, magnet_pos, mu, type) % Inputs are all in world frame, sens_or is a quaternion
    field_world_frame = mag_field(sens_pos-magnet_pos, mu);

    % Convert to sensor frame
    sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(sens_or); % with respect to world frame
    field_sensor_frame = sensor_rotation_matrix.' * field_world_frame;
    
    if type == "3D"
        output = field_sensor_frame;
    elseif type == "1D"
        output = field_sensor_frame(3);
    end
end

function output = sensor_output_from_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type)
    Rx = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];

    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];

    mu_world = Rx * Ry * [0;0;mu_norm];

    B_world = mag_field(sens_pos-magnet_pos, mu_world);

    % Convert to sensor frame
    sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(sens_or); % sensor to world frame
    B_sensor = sensor_rotation_matrix.' * B_world;

    if type == "3D"
        output = B_sensor;
    elseif type == "1D"
        output = B_sensor(3);
    end
end

function output = sensor_output_from_quaternion(sens_pos, sens_or, magnet_pos, mu_norm, w, x, y, z, type)
    q = [w;x;y;z];
    
    R = quaternionToMatrix(q); % magnet to world

    mu_world = R * [0;0;mu_norm];

    B_world = mag_field(sens_pos-magnet_pos, mu_world);

    % Convert to sensor frame
    sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(sens_or); % sensor to world frame
    B_sensor = sensor_rotation_matrix.' * B_world;

    if type == "3D"
        output = B_sensor;
    elseif type == "1D"
        output = B_sensor(3);
    end
end


% Jacobian numerical
% From B in sensor frame to magnet position in world frame
function J = jacobian_numerical(sens_pos, sens_or, magnet_pos, mu, delta, type)
    % ps: sensor position, pm: magnet position, mu: magnetic dipole moment
    % delta: increment
    % B(ps) = f(pm, ps, mu) Compute the Jacobian with respect to pm
    % J(ps) = d f(pm, ps) / d pm
    
    senser_out = sensor_output(sens_pos, sens_or, magnet_pos, mu, type);
    
    % Numerical derivative
    delta_x = magnet_pos + [delta;0;0];
    delta_y = magnet_pos + [0;delta;0];
    delta_z = magnet_pos + [0;0;delta];
    
    
    J = [(sensor_output(sens_pos, sens_or, delta_x, mu, type)-senser_out)/delta, ... 
         (sensor_output(sens_pos, sens_or, delta_y, mu, type)-senser_out)/delta, ... 
         (sensor_output(sens_pos, sens_or, delta_z, mu, type)-senser_out)/delta];
end

% Jacobian numerical
% From B in sensor frame to magnet position in world frame
function J = jacobian_numerical_sensor_reading_to_magnet_conf_euler_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, delta, type)
    % ps: sensor position, pm: magnet position, mu: magnetic dipole moment
    % delta: increment
    % B(ps) = f(pm, ps, mu) Compute the Jacobian with respect to pm
    % J(ps) = d f(pm, ps) / d pm
    
    senser_out = sensor_output_from_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type);
    
    % Numerical derivative
    delta_x = magnet_pos + [delta;0;0];
    delta_y = magnet_pos + [0;delta;0];
    delta_z = magnet_pos + [0;0;delta];

    delta_theta = theta + delta;
    delta_phi = phi + delta;
    
    
    J = [(sensor_output_from_angle(sens_pos, sens_or, delta_x, mu_norm, theta, phi, type)-senser_out)/delta, ... 
         (sensor_output_from_angle(sens_pos, sens_or, delta_y, mu_norm, theta, phi, type)-senser_out)/delta, ... 
         (sensor_output_from_angle(sens_pos, sens_or, delta_z, mu_norm, theta, phi, type)-senser_out)/delta, ...
         (sensor_output_from_angle(sens_pos, sens_or, magnet_pos, mu_norm, delta_theta, phi, type)-senser_out)/delta, ...
         (sensor_output_from_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, delta_phi, type)-senser_out)/delta];
end

% Jacobian numerical
% From B in sensor frame to magnet position in world frame
function J = jacobian_numerical_sensor_reading_to_magnet_conf_quaternion(sens_pos, sens_or, magnet_pos, mu_norm, w, x, y, z, delta, type)
    % ps: sensor position, pm: magnet position, mu: magnetic dipole moment
    % delta: increment
    % B(ps) = f(pm, ps, mu) Compute the Jacobian with respect to pm
    % J(ps) = d f(pm, ps) / d pm
    
    senser_out = sensor_output_from_quaternion(sens_pos, sens_or, magnet_pos, mu_norm, w, x, y, z, type);
    
    % Numerical derivative
    delta_x = magnet_pos + [delta;0;0];
    delta_y = magnet_pos + [0;delta;0];
    delta_z = magnet_pos + [0;0;delta];

    delta_w = w + delta;
    delta_x_quaternion = x + delta;
    delta_y_quaternion = y + delta;
    delta_z_quaternion = z + delta;
    
    
    J = [(sensor_output_from_quaternion(sens_pos, sens_or, delta_x, mu_norm, w, x, y, z, type)-senser_out)/delta, ... 
         (sensor_output_from_quaternion(sens_pos, sens_or, delta_y, mu_norm, w, x, y, z, type)-senser_out)/delta, ... 
         (sensor_output_from_quaternion(sens_pos, sens_or, delta_z, mu_norm, w, x, y, z, type)-senser_out)/delta, ...
         (sensor_output_from_quaternion(sens_pos, sens_or, magnet_pos, mu_norm, delta_w, x, y, z, type)-senser_out)/delta, ...
         (sensor_output_from_quaternion(sens_pos, sens_or, magnet_pos, mu_norm, w, delta_x_quaternion, y, z, type)-senser_out)/delta, ...
         (sensor_output_from_quaternion(sens_pos, sens_or, magnet_pos, mu_norm, w, x, delta_y_quaternion, z, type)-senser_out)/delta, ...
         (sensor_output_from_quaternion(sens_pos, sens_or, magnet_pos, mu_norm, w, x, y, delta_z_quaternion, type)-senser_out)/delta];
end

% Jacobian analytical
% Jacobian of B in world frame to magnet position in world frame
function J = jacobian_analytical_world_to_world(sens_pos, sens_or, magnet_pos, mu, type)
    mu0 = 4*pi*1e-7;
    r = sens_pos - magnet_pos;
    r_hat = r/norm(r);

    mu_hat = mu/norm(mu);

    unit_sens_or = sens_or/norm(sens_or);

    % Convert to sensor frame
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % with respect to world frame

    J = ((3*mu0*norm(mu))/(4*pi*norm(r)^4)) * ((eye(3)-5*r_hat*r_hat.')*(r_hat.'*mu_hat) + mu_hat*r_hat.' + r_hat*mu_hat.');
    
    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame
function J = jacobian_analytical_sensor_reading_to_world(sens_pos, sens_or, magnet_pos, mu, type)
    mu0 = 4*pi*1e-7;

    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r = sens_pos - magnet_pos;
    
    r = sensor_rotation_matrix.'*r; % convert from world frame to sensor frame
    r_hat = r/norm(r);
    
    mu = sensor_rotation_matrix.'*mu; % convert from world frame to sensor frame
    mu_hat = mu/norm(mu);

    J = -((3*mu0*norm(mu))/(4*pi*norm(r)^4)) * ((eye(3)-5*r_hat*r_hat.')*(r_hat.'*mu_hat) + mu_hat*r_hat.' + r_hat*mu_hat.');
    
    J = J * sensor_rotation_matrix.';
    
    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame, euler angle
function J = jacobian_B_sensor_to_magnet_conf(sens_pos, sens_or, magnet_pos, mu, theta, phi, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r = sens_pos - magnet_pos;
    
    r = sensor_rotation_matrix.'*r; % convert from world frame to sensor frame
    r_hat = r/norm(r);
    
    mu_sensor = sensor_rotation_matrix.'*mu; % convert from world frame to sensor frame
    mu_sensor_hat = mu_sensor/norm(mu_sensor);

    J_position = -((3*mu0*norm(mu_sensor))/(4*pi*norm(r)^4)) * ((eye(3)-5*(r_hat*r_hat.'))*(r_hat.'*mu_sensor_hat) + mu_sensor_hat*r_hat.' + r_hat*mu_sensor_hat.');
    
    J_position = J_position * sensor_rotation_matrix.';

    % Second part of the jacobian, from sensor reading to magnet
    % orientation

    % Extract rotation information
    mu_norm = norm(mu);
    dmu_dangle = mu_norm*[0 cos(phi); -cos(theta)*cos(phi) sin(theta)*sin(phi); -sin(theta)*cos(phi) -cos(theta)*sin(phi)];

    J_angle = (mu0/(4*pi))*(1/(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*sensor_rotation_matrix.'*dmu_dangle;

    J = [J_position, J_angle];

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end

    
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame, quaternion
function J = jacobian_B_sensor_to_magnet_conf_quaternion(sens_pos, sens_or, magnet_pos, mu, w, x, y ,z , type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r = sens_pos - magnet_pos;
    
    r = sensor_rotation_matrix.'*r; % convert from world frame to sensor frame
    r_hat = r/norm(r);
    
    mu_sensor = sensor_rotation_matrix.'*mu; % convert from world frame to sensor frame
    mu_sensor_hat = mu_sensor/norm(mu_sensor);

    J_position = -((3*mu0*norm(mu_sensor))/(4*pi*norm(r)^4)) * ((eye(3)-5*(r_hat*r_hat.'))*(r_hat.'*mu_sensor_hat) + mu_sensor_hat*r_hat.' + r_hat*mu_sensor_hat.');
    
    J_position = J_position * sensor_rotation_matrix.';

    % Second part of the jacobian, from sensor reading to magnet
    % orientation

    % Extract rotation information
    mu_norm = norm(mu);
    dmu_dquaternion = mu_norm * [2*y, 2*z, 2*w, 2*x; -2*x, -2*w, 2*z, 2*y; 0 -4*x, -4*y, 0];

    J_angle = (mu0/(4*pi))*(1/(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*sensor_rotation_matrix.'*dmu_dquaternion;

    J = [J_position, J_angle];

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end

    
end



%% Objective function
function obj = min_fun(sens_conf, magnet_conf, mu, type)
    
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    sens_conf = reshape(sens_conf, 7, []);
    
    % Extract position and orientation
    sens_pos = sens_conf(1:3,:);

    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Collect the reciprocal condition number and the min magnetic field strength for each magnet configuration
    reciprocal_condition_number_set = [];
    min_singular_value = [];
    min_magnetic_field_set = [];
    for magnet_num=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(:,magnet_num);
        
        J = [];
        for i=1:sens_num
            J = [J;jacobian_analytical_sensor_reading_to_world(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type)];
        end

        sigma = svd(J);

        num_dof = size(J,2);
        reciprocal_condition_number = sigma(num_dof)/sigma(1);
        
        % Construct the minimum sensor outputs for each sensor in the list
        sensor_outputs = [];
        for i=1:sens_num
            sensor_output_one_sensor = sensor_output(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type);
            % Put the min absolute sensor reading into the list, we care
            % about the magnitude
            sensor_outputs = [sensor_outputs, min(abs(sensor_output_one_sensor))]; 
        end

        reciprocal_condition_number_set = [reciprocal_condition_number_set; reciprocal_condition_number]; % Put the rcond for this magnet conf into the list
        min_magnetic_field_set = [min_magnetic_field_set; min(sensor_outputs)]; % Put the min sensor reading of all sensors in to the list
        min_singular_value = [min_singular_value; sigma(num_dof)];
    end

    % Minimize the negative of the min in the list -> maximize the min in
    % the list
    obj = [-min(reciprocal_condition_number_set)];
end

%% Evaluation function
function [obj_rcond_all_solution, obj_B_all_solution, min_rcond, min_B, min_singular_value] = evaluate(sens_conf, magnet_conf, mu, type)
    optimal_sol_num = size(sens_conf, 1);
    obj_rcond_all_solution = [];
    obj_B_all_solution = [];
    obj_min_singular_value_all_solution = [];
    

    for num=1:optimal_sol_num
        sens_conf_one = sens_conf(num,:);
        sens_num = sens_conf_one(end);
        sens_conf_one(:,end) = [];
        sens_conf_one = reshape(sens_conf_one, 7, []);

        % Extract position and orientation
        sens_pos = sens_conf_one(1:3,:);
    
        sens_or = sens_conf_one(4:7,:);
        magnitudes = vecnorm(sens_or);
        sens_or_unitary = sens_or ./ magnitudes;

        reciprocal_number_all_magnet_conf = [];
        min_sensor_reading_all_magnet_conf = [];
        min_singular_value_all_magnet_conf = [];
        for magnet_num=1:size(magnet_conf,2)

            magnet_pos = magnet_conf(:,magnet_num);
            J = [];
            for i=1:sens_num
                J = [J;jacobian_analytical_sensor_reading_to_world(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type)];
            end
        
            sigma = svd(J);
            reciprocal_condition_number = sigma(end)/sigma(1);
            
            % Construct the minimum sensor outputs for each sensor in the list
            sensor_outputs = [];
            for i=1:sens_num
                sensor_output_one_sensor = sensor_output(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type);
                % Put the min absolute sensor reading into the list, we care
                % about the magnitude
                sensor_outputs = [sensor_outputs, min(abs(sensor_output_one_sensor))]; 
            end
            
            reciprocal_number_all_magnet_conf = [reciprocal_number_all_magnet_conf;reciprocal_condition_number];
            min_sensor_reading_all_magnet_conf = [min_sensor_reading_all_magnet_conf; min(sensor_outputs)];
            min_singular_value_all_magnet_conf = [min_singular_value_all_magnet_conf; sigma(end)];
        end

        obj_rcond_all_solution = [obj_rcond_all_solution, reciprocal_number_all_magnet_conf];
        obj_B_all_solution = [obj_B_all_solution, min_sensor_reading_all_magnet_conf];
        obj_min_singular_value_all_solution = [obj_min_singular_value_all_solution, min_singular_value_all_magnet_conf];
    end
    min_rcond = min(obj_rcond_all_solution);
    min_B = min(obj_B_all_solution);
    min_singular_value = min(obj_min_singular_value_all_solution);
end

%% LSQ Object function
function obj = min_fun_lsq(sens_conf, magnet_conf, mu)
    
    sens_conf = reshape(sens_conf, 3, []);
    sens_num = 4;
    J = [];
    for i=1:sens_num
        J = [J;jacobian_analytical(sens_conf(:,i), magnet_conf, mu)];
    end

    sigma = svd(J);
    reciprocal_condition_number = sigma(end)/sigma(1);
    
    B = [];
    for i=1:sens_num
        r = sens_conf(:,sens_num) - magnet_conf;
        B = [B, norm(mag_field(r, mu))];
    end
    
    obj = [-reciprocal_condition_number]
end