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
mu_norm = norm(mu);

delta = 1e-7;
type = "1D";

unit = 'millimeter';
if strcmp(unit, 'meter') 
    scale = 1;
elseif strcmp(unit, 'decimeter')
    scale = 0.1;
elseif strcmp(unit, 'centimeter')
    scale = 0.01;
elseif strcmp(unit, 'millimeter')
    scale = 0.001;
end

%% Magnet configurations 1
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

%% Magnet configurations 2
% Workspace as a plane in m
LB = [-0.05, -0.05];
UB = [0.05, 0.05];

% Number of samples within the box will be num_samp^3
num_samp = 2;

C = cell(1, 2);
[C{:}] = ndgrid(linspace(0, 1, num_samp));
C = cellfun(@(a) a(:), C, 'Uni', 0);
combo = [C{:}];
conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
    repmat(LB, [size(combo, 1), 1]);
conf = conf.';

magnet_conf = conf;

z_values = 0.2 * ones(1, size(magnet_conf,2));
magnet_conf = [magnet_conf;z_values];

%% Magnet configurations with orientation
% Workspace as a plane in m
x = [0];
y = [0];
z = [0.2/scale];
theta = [0];
phi = [0];

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


%% Genetic algorithm
lb = [-0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      -0.2/scale -0.2/scale 0.0 -1 -1 -1 -1 ...
      5];
ub = [0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      0.2/scale 0.2/scale 0.0 1 1 1 1 ...
      5];
options = optimoptions(@gamultiobj,'Display','iter', 'MaxStallGenerations', 30000, 'MaxGenerations', 30000, 'PopulationSize', 30000);

fun = @(x) min_fun_orientation(x, magnet_conf, mu_norm, type);
[sol, fval, exitflag, output] = gamultiobj(fun, length(lb), [], [], [], [], lb, ub, [], length(lb), options);

save('results_5_one_axis_centimeter_max_rcond_4_30000gen_30000pop')
%%
% Evaluate
load('results_5_one_axis_centimeter_max_rcond_2_10000gen_10000pop.mat')
sens_conf = [sol];
[obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);
sol
obj_rcond
min_rcond

%%
mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

max(mean_rcond)
max(min_rcond)

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
% magnet orientation with respect to world frame, quaternion
function J = jacobian_analytical_sensor_reading_to_magnet_conf_quaternion(sens_pos, sens_or, magnet_pos, mu, w, x, y ,z , type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r  = sens_pos - magnet_pos;
    
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
function obj = min_fun_orientation(sens_conf, magnet_conf, mu_norm, type)
    
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
    min_singular_value_set = [];
%     min_magnetic_field_set = [];

    for magnet_num=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        
        J = [];
        for i=1:sens_num
            J = [J;jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle(sens_pos(:,i), sens_or_unitary(:,i), ...
                magnet_pos, mu_norm, theta, phi, type)];
        end

        sigma = svd(J);

        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);

        reciprocal_condition_number_set = [reciprocal_condition_number_set; reciprocal_condition_number]; % Put the rcond for this magnet conf into the list
        min_singular_value_set = [min_singular_value_set; sigma(num_dof)];
    end

%     % Minimize the negative of the min in the list -> maximize the min in
%     % the list
%     obj = [-min(reciprocal_condition_number_set)];
    
%     % Maximize all the reciprocal condition number
%     obj = [];
%     for i = 1:size(reciprocal_condition_number_set,1)
%         obj = [obj, -reciprocal_condition_number_set(i)];
%     end

    obj = [-min(reciprocal_condition_number_set)];
end


%% Evaluation function
function [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type)
    optimal_sol_num = size(sens_conf, 1);
    obj_rcond = [];
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

        reciprocal_number_one_magnet_conf = [];
        for magnet_num=1:size(magnet_conf,2)

            magnet_pos = magnet_conf(1:3,magnet_num);
            theta = magnet_conf(4, magnet_num);
            phi = magnet_conf(5, magnet_num);
            
            J = [];
            for i=1:sens_num
                J = [J;jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle(sens_pos(:,i), sens_or_unitary(:,i), ...
                    magnet_pos, mu_norm, theta, phi, type)];
            end
    
            sigma = svd(J);
    
            num_dof = 5;
            reciprocal_condition_number = sigma(num_dof)/sigma(1);
            
            reciprocal_number_one_magnet_conf = [reciprocal_number_one_magnet_conf;reciprocal_condition_number];
        end

        obj_rcond = [obj_rcond, reciprocal_number_one_magnet_conf];
    end
    min_rcond = min(obj_rcond);
end