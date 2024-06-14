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
sph_dip = 13.2*(4/3)*pi*(norm(sph_dia)/2)^3/mu0; % spheres magnetic dipole
mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole
mu_norm = norm(mu);

B_r = 13.2;
Volumn = (4/3)*pi*(norm(sph_dia)/2)^3; % volumn in the unit of (specified unit)^3

delta = 1e-7;
type = "1D";

%% Magnet configurations with orientation
% Workspace as a plane in m
x = [0];
y = [0];
z = [0.2/scale];
theta = [0];
phi = [0];
psi = [0];

% Create grids for each dimension
[X, Y, Z, Theta, Phi] = ndgrid(x, y, z, theta, phi);

% Reshape the grids into column vectors
X = reshape(X, [], 1);
Y = reshape(Y, [], 1);
Z = reshape(Z, [], 1);
Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);

% Combine all dimensions into a single matrix
magnet_conf = [X, Y, Z, Theta, Phi].';

magnet_conf = [0 0 0.2/scale 0 0 0;... 
               0 0 0.2/scale deg2rad(90) deg2rad(0) 0; ...
               0 0 0.2/scale deg2rad(0) deg2rad(90) 0].';


%% Genetic algorithm
lb = [-0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      12];
ub = [0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      12];
options = optimoptions(@gamultiobj,'Display','iter', 'MaxStallGenerations', 2000, 'MaxGenerations', 2000, 'PopulationSize', 20000);

fun = @(x) min_fun_orientation(x, magnet_conf, B_r, Volumn, type);
constraint_fun = @(x) jacobian_constraint(x, magnet_conf, mu_norm, type);
[sol, fval, exitflag, output] = gamultiobj(fun, length(lb), [], [], [], [], lb, ub, [], length(lb), options);

save('results_12_1_axis_multiobj_2000gen_20000pop_scaled_unit')
%%
% Evaluate
load('results_12_1_axis_multiobj_2000gen_20000pop_scaled_unit.mat')
sens_conf = [sol];
[obj_rcond, min_rcond, obj_svd, min_svd] = evaluate_with_orientation(sens_conf, magnet_conf, B_r, Volumn, type);
sol
obj_rcond
min_rcond
obj_svd
min_svd
%%
scatter(obj_rcond, obj_svd)
%% Jacobian test
sens_conf = [0.2 0.2 0 1 0 0 0 ...
             0.2 -0.2 0 1 0 0 0 ...
             -0.2 0.2 0 1 0 0 0 ...
             -0.2 -0.2 0 1 0 0 0 ...
             4];
magnet_conf = [0;0;0.2;deg2rad(0);deg2rad(0);deg2rad(0)];

[obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

sens_num = sens_conf(end);
sens_conf(end) = [];
sens_conf = reshape(sens_conf, 7, []);

% Extract position and orientation
sens_pos = sens_conf(1:3,:);

sens_or = sens_conf(4:7,:);
magnitudes = vecnorm(sens_or);
sens_or_unitary = sens_or ./ magnitudes;

% Collect the reciprocal condition number for each magnet configuration
for magnet_num=1:size(magnet_conf,2)
    magnet_pos = magnet_conf(1:3,magnet_num);
    theta = magnet_conf(4, magnet_num);
    phi = magnet_conf(5, magnet_num);
    psi = magnet_conf(6, magnet_num);

    J = [];
    for i=1:sens_num
        J = [J;jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_XYX(sens_pos(:,i), sens_or_unitary(:,i), ...
            magnet_pos, mu_norm, theta, phi, psi, type)];
    end

    sigma = svd(J);

    num_dof = 5;
    reciprocal_condition_number = sigma(num_dof)/sigma(1);
end

%% Jacobian test - scaled unit
scale = 0.15;

sens_conf = [0.5/scale 0.5/scale 0 1 0 0 0 ...
             0.5/scale -0.5/scale 0 1 0 0 0 ...
             -0.5/scale 0.5/scale 0 1 0 0 0 ...
             -0.5/scale -0.5/scale 0 1 0 0 0 ...
             4];
magnet_conf = [0;0;0.2/scale;deg2rad(0);deg2rad(0);deg2rad(0)];

Volumn = Volumn / (scale^3);

sens_num = sens_conf(end);
sens_conf(end) = [];
sens_conf = reshape(sens_conf, 7, []);

% Extract position and orientation
sens_pos = sens_conf(1:3,:);

sens_or = sens_conf(4:7,:);
magnitudes = vecnorm(sens_or);
sens_or_unitary = sens_or ./ magnitudes;

% Collect the reciprocal condition number for each magnet configuration
for magnet_num=1:size(magnet_conf,2)
    magnet_pos = magnet_conf(1:3,magnet_num);
    theta = magnet_conf(4, magnet_num);
    phi = magnet_conf(5, magnet_num);
    psi = magnet_conf(6, magnet_num);

    J_scaled = [];
    for i=1:sens_num
        J_scaled = [J_scaled;J_analytical_sensor_B_to_magnet_conf_euler_angle_XYX_scaled_uni(sens_pos(:,i), sens_or_unitary(:,i), ...
            magnet_pos, B_r, Volumn, theta, phi, psi, type)];
    end

    sigma_scaled = svd(J_scaled);

    num_dof = 5;
    reciprocal_condition_number_scaled = sigma_scaled(num_dof)/sigma_scaled(1);
end

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

function output = sensor_output_from_euler_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type)
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
    sensor_rotation_matrix = quaternionToMatrix(sens_or); % with respect to world frame
    B_world = sensor_rotation_matrix.' * B_world;

    if type == "3D"
        output = B_world;
    elseif type == "1D"
        output = B_world(3);
    end
end

function output = sensor_output_from_quaternion(sens_pos, sens_or, magnet_pos, mu_norm, w, x, y, z, type)
    q = [w;x;y;z];
    
    R = quaternionToMatrix(q);

    mu_world = R * [0;0;mu_norm];

    B_world = mag_field(sens_pos-magnet_pos, mu_world);

    % Convert to sensor frame
    sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(sens_or); % with respect to world frame
    B_world = sensor_rotation_matrix.' * B_world;

    if type == "3D"
        output = B_world;
    elseif type == "1D"
        output = B_world(3);
    end
end

function output = sensor_output_from_spherical_ang(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type)
    mu_world = mu_norm*[cos(theta)*cos(phi);cos(theta)*sin(phi);sin(theta)];

    B_world = mag_field(sens_pos-magnet_pos, mu_world);

    % Convert to sensor frame
    sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(sens_or); % with respect to world frame
    B_world = sensor_rotation_matrix.' * B_world;

    if type == "3D"
        output = B_world;
    elseif type == "1D"
        output = B_world(3);
    end
end


%% Jacobian numerical
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
    
    senser_out = sensor_output_from_euler_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type);
    
    % Numerical derivative
    delta_x = magnet_pos + [delta;0;0];
    delta_y = magnet_pos + [0;delta;0];
    delta_z = magnet_pos + [0;0;delta];

    delta_theta = theta + delta;
    delta_phi = phi + delta;
    
    
    J = [(sensor_output_from_euler_angle(sens_pos, sens_or, delta_x, mu_norm, theta, phi, type)-senser_out)/delta, ... 
         (sensor_output_from_euler_angle(sens_pos, sens_or, delta_y, mu_norm, theta, phi, type)-senser_out)/delta, ... 
         (sensor_output_from_euler_angle(sens_pos, sens_or, delta_z, mu_norm, theta, phi, type)-senser_out)/delta, ...
         (sensor_output_from_euler_angle(sens_pos, sens_or, magnet_pos, mu_norm, delta_theta, phi, type)-senser_out)/delta, ...
         (sensor_output_from_euler_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, delta_phi, type)-senser_out)/delta];
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

% Jacobian numerical
% From B in sensor frame to magnet position in world frame
% And magnet orientation expressed in spherical coordinate
function J = jacobian_numerical_sensor_reading_to_magnet_conf_spherical_ang(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, delta, type)
    % ps: sensor position, pm: magnet position, mu: magnetic dipole moment
    % delta: increment
    % B(ps) = f(pm, ps, mu) Compute the Jacobian with respect to pm
    % J(ps) = d f(pm, ps) / d pm
    
    senser_out = sensor_output_from_spherical_ang(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type);
    
    % Numerical derivative
    delta_x = magnet_pos + [delta;0;0];
    delta_y = magnet_pos + [0;delta;0];
    delta_z = magnet_pos + [0;0;delta];

    delta_theta = theta + delta;
    delta_phi = phi + delta;
    
    
    J = [(sensor_output_from_spherical_ang(sens_pos, sens_or, delta_x, mu_norm, theta, phi, type)-senser_out)/delta, ... 
         (sensor_output_from_spherical_ang(sens_pos, sens_or, delta_y, mu_norm, theta, phi, type)-senser_out)/delta, ... 
         (sensor_output_from_spherical_ang(sens_pos, sens_or, delta_z, mu_norm, theta, phi, type)-senser_out)/delta, ...
         (sensor_output_from_spherical_ang(sens_pos, sens_or, magnet_pos, mu_norm, delta_theta, phi, type)-senser_out)/delta, ...
         (sensor_output_from_spherical_ang(sens_pos, sens_or, magnet_pos, mu_norm, theta, delta_phi, type)-senser_out)/delta];
end

% From B in sensor frame to magnet position and mu in world frame
function J = jacobian_numerical_sensor_reading_to_magnet_conf_mu(sens_pos, sens_or, magnet_pos, mu, delta, type)
    % ps: sensor position, pm: magnet position, mu: magnetic dipole moment
    % delta: increment
    % B(ps) = f(pm, ps, mu) Compute the Jacobian with respect to pm
    % J(ps) = d f(pm, ps) / d pm
    
    senser_out = sensor_output(sens_pos, sens_or, magnet_pos, mu, type);
    
    % Numerical derivative
    delta_x = magnet_pos + [delta;0;0];
    delta_y = magnet_pos + [0;delta;0];
    delta_z = magnet_pos + [0;0;delta];

    delta_mux = mu + [delta;0;0];
    delta_muy = mu + [0;delta;0];
    delta_muz = mu + [0;0;delta];
    
    
    J = [(sensor_output(sens_pos, sens_or, delta_x, mu, type)-senser_out)/delta, ... 
         (sensor_output(sens_pos, sens_or, delta_y, mu, type)-senser_out)/delta, ... 
         (sensor_output(sens_pos, sens_or, delta_z, mu, type)-senser_out)/delta, ...
         (sensor_output(sens_pos, sens_or, magnet_pos, delta_mux, type)-senser_out)/delta, ...
         (sensor_output(sens_pos, sens_or, magnet_pos, delta_muy, type)-senser_out)/delta, ...
         (sensor_output(sens_pos, sens_or, magnet_pos, delta_muz, type)-senser_out)/delta];
end

%% Jacobian analytical
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
function J = jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    Rx = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];
    
    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];

    mu_world = Rx * Ry * [0;0;mu_norm];
    mu_world_hat = mu_world/norm(mu_world);

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r_world = sens_pos - magnet_pos;
    r_world_hat = r_world/norm(r_world);
    
    J_position = -((3*mu0*norm(mu_world))/(4*pi*norm(r_world)^4)) * ((eye(3)-5*(r_world_hat*r_world_hat.'))*(r_world_hat.'*mu_world_hat) ...
                + mu_world_hat*r_world_hat.' + r_world_hat*mu_world_hat.');

    % Second part of the jacobian, from sensor reading to magnet
    % orientation

    % Extract rotation information
    dmu_dangle = mu_norm*[0 cos(phi); -cos(theta)*cos(phi) sin(theta)*sin(phi); -sin(theta)*cos(phi) -cos(theta)*sin(phi)];

    J_angle = (mu0/(4*pi))*(1/(norm(r_world)^3))*(3*(r_world_hat*r_world_hat.')-eye(3))*dmu_dangle;

    J = [J_position, J_angle];

    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame, euler angle
function J = jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_XYX(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, psi, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    Rx = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];
    
    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];

    Rx_2 = [1 0 0;                  % rotate about x-axis
    0 cos(psi) -sin(psi);
    0 sin(psi) cos(psi)];

    mu_world = Rx * Ry * Rx_2 * [0;0;mu_norm];
    mu_world_hat = mu_world/norm(mu_world);

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r_world = sens_pos - magnet_pos;
    r_world_hat = r_world/norm(r_world);
    
    J_position = -((3*mu0*norm(mu_world))/(4*pi*norm(r_world)^4)) * ((eye(3)-5*(r_world_hat*r_world_hat.'))*(r_world_hat.'*mu_world_hat) ...
                + mu_world_hat*r_world_hat.' + r_world_hat*mu_world_hat.');

    % Second part of the jacobian, from sensor reading to magnet
    % orientation

    % Extract rotation information
    dmu_dangle = mu_norm*[0 cos(phi)*cos(psi) -sin(phi)*sin(psi); 
                          sin(theta)*sin(psi)-cos(phi)*cos(psi)*cos(theta) sin(phi)*cos(psi)*sin(theta) -cos(theta)*cos(psi)+cos(phi)*sin(psi)*sin(theta); 
                          -cos(phi)*cos(psi)*sin(theta)-sin(psi)*cos(theta) -sin(phi)*cos(psi)*cos(theta) -cos(phi)*sin(psi)*cos(theta)-cos(psi)*sin(theta)];

    J_angle = (mu0/(4*pi))*(1/(norm(r_world)^3))*(3*(r_world_hat*r_world_hat.')-eye(3))*dmu_dangle;

    J = [J_position, J_angle];

    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
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

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame, quaternion
function J = jacobian_analytical_sensor_reading_to_magnet_conf_quaternion(sens_pos, sens_or, magnet_pos, mu, w, x, y ,z , type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
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
    dmu_dquaternion = mu_norm * [2*y, 2*z, 2*w, 2*x; -2*x, -2*w, 2*z, 2*y; 0, -4*x, -4*y, 0];

    J_angle = (mu0/(4*pi))*(1/(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*sensor_rotation_matrix.'*dmu_dquaternion;

    J = [J_position, J_angle];

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end 
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame, quaternion
function J = jacobian_analytical_sensor_reading_to_magnet_conf_mu(sens_pos, sens_or, magnet_pos, mu, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r_world = sens_pos - magnet_pos;
    r_world_hat = r_world/norm(r_world);
    
    mu_world = mu;
    mu_world_hat = mu_world/norm(mu_world);

    J_position = -((3*mu0*norm(mu_world))/(4*pi*norm(r_world)^4)) * ((eye(3)-5*(r_world_hat*r_world_hat.'))*(r_world_hat.'*mu_world_hat) ...
                    + mu_world_hat*r_world_hat.' + r_world_hat*mu_world_hat.');

    % Second part of the jacobian, from sensor reading to magnet
    % orientation, embedded in mu

    J_mu = (mu0/(4*pi)) * (1/(norm(r_world)^3)) * (3*r_world_hat*r_world_hat.'-eye(3));

    J = [J_position, J_mu];
    
    % Convert to sensor frame
    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end 
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame, quaternion
function J = jacobian_analytical_sensor_reading_to_magnet_conf_spherical_ang(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r_world = sens_pos - magnet_pos;
    r_world_hat = r_world/norm(r_world);
    
    mu_world = mu_norm*[cos(theta)*cos(phi); cos(theta)*sin(phi); sin(theta)];
    mu_world_hat = mu_world/norm(mu_world);

    J_position = -((3*mu0*mu_norm)/(4*pi*norm(r_world)^4)) * ((eye(3)-5*(r_world_hat*r_world_hat.'))*(r_world_hat.'*mu_world_hat) ...
                    + mu_world_hat*r_world_hat.' + r_world_hat*mu_world_hat.');

    % Second part of the jacobian, from sensor reading to magnet
    % orientation, spherical coordinate

    J_spherical_ang = (mu0/(4*pi)) * (1/(norm(r_world)^3)) * (3*r_world_hat*r_world_hat.'-eye(3)) * ...
                        mu_norm*[-sin(theta)*cos(phi) -cos(theta)*sin(phi); -sin(theta)*sin(phi) cos(theta)*cos(phi); cos(theta) 0];

    J = [J_position, J_spherical_ang];
    
    % Convert to sensor frame
    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end 
end

%% Objective function
function obj = min_fun_orientation(sens_conf, magnet_conf, B_r, Volumn, type)
    
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    sens_conf = reshape(sens_conf, 7, []);
    
    % Extract position and orientation
    sens_pos = sens_conf(1:3,:);

    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Collect the reciprocal condition number for each magnet configuration
    rcond_set = [];
    min_singular_value_set = [];
    min_singular_value_Jp_set = [];
    min_singular_value_Jo_set = [];

    for magnet_num=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
        
        J = [];
        for i=1:sens_num
            J = [J;J_analytical_sensor_B_to_magnet_conf_euler_angle_XYX_scaled_uni(sens_pos(:,i), sens_or_unitary(:,i), ...
                magnet_pos, B_r, Volumn, theta, phi, psi, type)];
        end

        sigma = svd(J);

        num_dof = 5;

        min_singular_value_set = [min_singular_value_set;sigma(num_dof)];
        rcond_set = [rcond_set; sigma(num_dof)/sigma(1)];

        Jp = J(:,1:3);
        Jo = J(:,4:6);
        sigma_Jp = svd(Jp);
        sigma_Jo = svd(Jo);

        min_singular_value_Jp_set = [min_singular_value_Jp_set; sigma_Jp(3)];
        min_singular_value_Jo_set = [min_singular_value_Jo_set; sigma_Jo(2)];
    end

    % Minimize the negative of the min in the list -> maximize the min in
    % the list
    obj = [-min(min_singular_value_Jp_set),-min(min_singular_value_Jo_set)];
    
%     % Maximize all the reciprocal condition number
%     obj = [];
%     for i = 1:size(min_singular_value_set,1)
%         obj = [obj, -min_singular_value_Jp_set(i), -min_singular_value_Jo_set(i)];
%     end

end

%% Jacobian constraint on Jp and Jo respectively
function [c, ceq] = jacobian_constraint(sens_conf, magnet_conf, mu_norm, type)
    % No equality constraint
    ceq = [];
    
    % Extract sensor information: # of sensors, sensor positions,
    % orientations
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    sens_conf = reshape(sens_conf, 7, []);
    
    % Extract position and orientation
    sens_pos = sens_conf(1:3,:);

    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;

    % Collect the reciprocal condition number for each magnet configuration
    % the rcond is calculated for the position Jacobian and orientation
    % Jacobian separately
    reciprocal_condition_number_Jp = [];
    reciprocal_condition_number_Jo = [];
    
    % For each magnet configuration, construct the Jp and Jo, get their
    % corresponding rcond
    for magnet_num=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
        
        % Concatenate the Jacobians for all the sensors
        J = [];
        for i=1:sens_num
            J = [J;jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_XYX(sens_pos(:,i), sens_or_unitary(:,i), ...
                magnet_pos, mu_norm, theta, phi, psi, type)];
        end
        
        % First 3 columns are Jp, the last 3 columns are Jo
        Jp = J(:,1:3);
        Jo = J(:,4:6);

        num_dof_p = 3;
        num_dof_o = 2;

        sigma_p = svd(Jp);
        sigma_o = svd(Jo);

        rcond_p = sigma_p(num_dof_p)/sigma_o(1);
        rcond_o = sigma_o(num_dof_o)/sigma_o(1);

        % Put the rconds into the list
        reciprocal_condition_number_Jp = [reciprocal_condition_number_Jp; rcond_p];
        reciprocal_condition_number_Jo = [reciprocal_condition_number_Jo; rcond_o];
    end

    % Extract the minimum rcond for Jp and Jo for all magnet configs
    min_rcond_p = min(reciprocal_condition_number_Jp);
    min_rcond_o = min(reciprocal_condition_number_Jo);
    
    % Construct the nonlinear constraints c(1) <= 0, c(2) <= 0
    % It becomes: min_rcond_p >= 0.5, min_rcond_o >= 0.5
    c(1) = -min_rcond_p + 0.5;
    c(2) = -min_rcond_o + 0.5;


end


%% Evaluation function
function [obj_rcond, min_rcond, obj_svd, min_svd] = evaluate_with_orientation(sens_conf, magnet_conf, B_r, Volumn, type)
    optimal_sol_num = size(sens_conf, 1);
    obj_rcond = [];
    obj_svd = [];
    obj_B = [];
    

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
        
        min_svd_one_magnet_conf = [];
        reciprocal_number_one_magnet_conf = [];
        for magnet_num=1:size(magnet_conf,2)

            magnet_pos = magnet_conf(1:3,magnet_num);
            theta = magnet_conf(4, magnet_num);
            phi = magnet_conf(5, magnet_num);
            psi = magnet_conf(6, magnet_num);
            
            J = [];
            for i=1:sens_num
                J = [J;J_analytical_sensor_B_to_magnet_conf_euler_angle_XYX_scaled_uni(sens_pos(:,i), sens_or_unitary(:,i), ...
                    magnet_pos, B_r, Volumn, theta, phi, psi, type)];
            end
    
            sigma = svd(J);
    
            num_dof = 5;
            reciprocal_condition_number = sigma(num_dof)/sigma(1);
            
            reciprocal_number_one_magnet_conf = [reciprocal_number_one_magnet_conf;reciprocal_condition_number];
            min_svd_one_magnet_conf = [min_svd_one_magnet_conf; sigma(num_dof)];
        end

        obj_rcond = [obj_rcond, reciprocal_number_one_magnet_conf];
        obj_svd = [obj_svd, min_svd_one_magnet_conf];
    end
    min_rcond = min(obj_rcond);
    min_svd = min(obj_svd);
end