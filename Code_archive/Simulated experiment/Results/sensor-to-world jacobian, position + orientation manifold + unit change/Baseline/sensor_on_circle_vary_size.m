%%
close all
clear
clc

%%
% Spheres
unit = 'decimeter';
if strcmp(unit, 'meter') 
    scale =1;
elseif strcmp(unit, 'decimeter')
    scale = 0.1;
elseif strcmp(unit, 'centimeter')
    scale = 0.01;
elseif strcmp(unit, 'millimeter')
    scale = 0.001;
end

mu0 =  4*pi*1e-7; % air permeability
sph_dia = 3.175e-3*[0; 0; 1] / scale; % sphere diameter 3.175e-3 m (1/8 inch)
                                      % converted to corresponding unit
                                      % using the scale factor

B_r = 1.32; % Residual flux density (T)
Volumn = (4/3)*pi*(norm(sph_dia)/2)^3; % volumn of the sphere dipole
                                       % in the unit of (specified unit)^3
                                       
sph_dip = B_r*Volumn/mu0; % spheres magnetic dipole mu_norm = B_r*V/mu0

mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole moment vector
mu_norm = norm(mu);

delta = 1e-7;
type = "3D";

%% Magnet configurations with orientation
% Workspace as a plane in unit
x = [-0.05, 0, 0.05]/scale;
y = [-0.05, 0, 0.05]/scale;
z = [0.05, 0.10, 0.15] /scale;
theta = [-pi/2, -pi/4, 0, pi/4, pi/2];
phi = [-pi/2, -pi/4, 0, pi/4, pi/2];
psi = [0];

% Create grids for each dimension
[X, Y, Z, Theta, Phi, Psi] = ndgrid(x, y, z, theta, phi, psi);

% Reshape the grids into column vectors
X = reshape(X, [], 1);
Y = reshape(Y, [], 1);
Z = reshape(Z, [], 1);
Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);
Psi = reshape(Psi, [], 1);

% Combine all dimensions into a single matrix
magnet_conf = [X, Y, Z, Theta, Phi, Psi].';

%% 
% Define the base radius increment
l = 0.002; % You can change this to any value of the base radius

% Number of circles
n = 10;

% Number of sets
num_sets = 1000;

% Initialize cell array to hold all sets of coordinates
sensor_config_collection = cell(num_sets, 1);

figure;
hold on;

% Loop through each set
for set_idx = 1:num_sets
    % Base radius for the current set
    base_radius = set_idx * l;
    
    % Initialize arrays to hold coordinates for the current set
    coordinates_set = cell(10, 1);
    
    % Loop through each pair of circles in the set
    for k = 1:5
        % Calculate the radii for the current pair
        radius1 = base_radius * k;
        radius2 = base_radius * (k + 0.5);
        
        % Angles for the points on the first circle
        angles1 = linspace(0, 2*pi, n+1);
        angles1(end) = []; % Remove the last element to get 10 points
        
        % Calculate coordinates for the first circle
        x1 = radius1 * cos(angles1);
        y1 = radius1 * sin(angles1);
        
        % Calculate midpoints angles for the second circle
        angles2 = angles1 + (pi/n);
        
        % Calculate coordinates for the second circle
        x2 = radius2 * cos(angles2);
        y2 = radius2 * sin(angles2);
        
        % Store the coordinates for the current pair of circles
        coordinates_set{2*k-1} = [x1', y1'];
        coordinates_set{2*k} = [x2', y2'];
        
        % Plot the circles and points (optional, for visual confirmation)
        if set_idx == 1 || set_idx == length(sensor_config_collection) % Only plot for the first set to avoid clutter
            theta = linspace(0, 2*pi, 100);
            circle1_x = radius1 * cos(theta);
            circle1_y = radius1 * sin(theta);
            circle2_x = radius2 * cos(theta);
            circle2_y = radius2 * sin(theta);
            
            plot(circle1_x, circle1_y, 'b'); % Plot the first circle
            plot(x1, y1, 'ro'); % Plot the points on the first circle
            plot(circle2_x, circle2_y, 'g'); % Plot the second circle
            plot(x2, y2, 'mo'); % Plot the points on the second circle
        end
    end
    
    % Store the coordinates set in the sensor_config_collection cell array
    sensor_config_collection{set_idx} = coordinates_set;
end

axis equal;
title('Equally Distant Points on Multiple Circles with Increasing Base Radii');
xlabel('x');
ylabel('y');
grid on;

% Display the coordinates for all sets
for set_idx = 1:num_sets
    disp(['Coordinates for set ', num2str(set_idx), ':']);
    for circle_idx = 1:10
        disp(['Coordinates for circle ', num2str(circle_idx), ':']);
        disp(sensor_config_collection{set_idx}{circle_idx});
    end
end

%% Iterate through the sensor config collection
rcond_all_sens_conf = [];
min_svd_all_sens_conf = [];
max_svd_all_sens_conf = [];

for i = 1:length(sensor_config_collection)
    sensor_config_one = sensor_config_collection{i};
    sens_num = 100;
    sens_pos = [];
    for circle_idx = 1:10
        sens_pos = [sens_pos;sensor_config_one{circle_idx}];
    end
    sens_pos = sens_pos.';

    % Add the line of zero as the z coordinate to the position of sensors
    z_coordinate = zeros(1, size(sens_pos, 2));
    sens_pos = [sens_pos; z_coordinate];
    sens_pos = sens_pos/scale;

    % Default orientation
    default_or = [1;0;0;0];

    % Orientation for all sensors
    sens_or_unitary = repmat(default_or, 1, sens_num);

    rcond_one_sens_conf = [];
    min_svd_one_sens_conf = [];
    max_svd_one_sens_conf = [];

    % Collect the reciprocal condition number for each magnet configuration
    for magnet_num=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
    
        Rx_1 = [1 0 0;                  % rotate about x-axis
        0 cos(theta) -sin(theta);
        0 sin(theta) cos(theta)];
        
        Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
        0 1 0;
        -sin(phi) 0 cos(phi)];
    
        Rx_2 = [1 0 0;                  % rotate about x-axis
        0 cos(psi) -sin(psi);
        0 sin(psi) cos(psi)];
    
        R_star = Rx_1 * Ry * Rx_2;
    
        J_scaled = [];
        for j=1:sens_num
            J_scaled = [J_scaled;J_analytical_sensor_B_to_world_pR(sens_pos(:,j), sens_or_unitary(:,j), ...
                magnet_pos, B_r, Volumn, R_star, type)];
        end
    
        sigma_scaled = svd(J_scaled);
    
        num_dof = 5;
        reciprocal_condition_number_scaled = sigma_scaled(num_dof)/sigma_scaled(1);
        
        rcond_one_sens_conf = [rcond_one_sens_conf; reciprocal_condition_number_scaled];
        min_svd_one_sens_conf = [min_svd_one_sens_conf; sigma_scaled(num_dof)];
        max_svd_one_sens_conf = [max_svd_one_sens_conf; sigma_scaled(1)];
    end
    
    % Put the min value in the workspace into the list
    rcond_all_sens_conf = [rcond_all_sens_conf; min(rcond_one_sens_conf)];
    min_svd_all_sens_conf = [min_svd_all_sens_conf; min(min_svd_one_sens_conf)];
    max_svd_all_sens_conf = [max_svd_all_sens_conf; min(max_svd_one_sens_conf)];
end

%% Plots
inner_radius = linspace(0.002, 2, 1000);
outer_radius = linspace(0.011, 11, 1000);
figure(1)
plot(outer_radius, rcond_all_sens_conf, "linewidth", 2);
% ylim([0 1])
xlim([0 11])
xlabel("radius of the outter circle (m)")
ylabel("Minimum reciprocal condition number in the workspace")
title("Minimum reciprocal condition number in the workspace vs side length of the square")
grid on
figure(2)
plot(inner_radius, rcond_all_sens_conf, "linewidth", 2);
% ylim([0 1])
xlim([0 2])
xlabel("radius of the inner circle (m)")
ylabel("Minimum reciprocal condition number in the workspace")
title("Minimum reciprocal condition number in the workspace vs side length of the square")
grid on

figure(3)
plot(outer_radius, min_svd_all_sens_conf, "linewidth", 2);
% ylim([0 2.5e-7])
xlim([0 11])
xlabel("radius of the outter circle (m)")
ylabel("Minimum of the \sigma_{min} in the workspace")
title("Minimum of the \sigma_{min} in the workspace vs side length of the square")
grid on
figure(4)
plot(inner_radius, min_svd_all_sens_conf, "linewidth", 2);
% ylim([0 2.5e-7])
xlim([0 2])
xlabel("radius of the inner circle (m)")
ylabel("Minimum of the \sigma_{min} in the workspace")
title("Minimum of the \sigma_{min} in the workspace vs side length of the square")
grid on

% figure(3)
% plot(2*side_length, max_svd_all_sens_conf, "linewidth", 2);
% % ylim([0 2.5e-7])
% xlabel("Side length of the square (m)")
% ylabel("Minimum of the \sigma_{max} in the workspace")
% title("Minimum of the \sigma_{max} in the workspace vs side length of the square")
% grid on


%% Jacobian test - scaled unit
scale = 0.1;

sens_conf = [0.05654/scale  0.05654/scale 0 1 0 0 0 ...
             0.05654/scale -0.05654/scale 0 1 0 0 0 ...
             -0.05654/scale 0.05654/scale 0 1 0 0 0 ...
             -0.05654/scale -0.05654/scale 0 1 0 0 0 ...
             ];

for magnet_num = 1:length(magnet_conf)
    sens_num = 4;
    sens_conf = reshape(sens_conf, 7, []);
    
    % Extract position and orientation
    sens_pos = sens_conf(1:3,:);
    
    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;

    rcond = [];
    min_svd = [];
    
    % Collect the reciprocal condition number for each magnet configuration
    for magnet_num=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
    
        Rx_1 = [1 0 0;                  % rotate about x-axis
        0 cos(theta) -sin(theta);
        0 sin(theta) cos(theta)];
        
        Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
        0 1 0;
        -sin(phi) 0 cos(phi)];
    
        Rx_2 = [1 0 0;                  % rotate about x-axis
        0 cos(psi) -sin(psi);
        0 sin(psi) cos(psi)];
    
        R_star = Rx_1 * Ry * Rx_2;
    
        J_scaled = [];
        for i=1:sens_num
            J_scaled = [J_scaled;J_analytical_sensor_B_to_world_pR(sens_pos(:,i), sens_or_unitary(:,i), ...
                magnet_pos, B_r, Volumn, R_star, type)];
        end
    
        sigma_scaled = svd(J_scaled);
    
        num_dof = 5;
        reciprocal_condition_number_scaled = sigma_scaled(num_dof)/sigma_scaled(1);
        rcond = [rcond; reciprocal_condition_number_scaled];
        min_svd = [min_svd; sigma_scaled(num_dof)];
    end
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

% Convert a vector to a skew symmetric matrix
function v_hat = skew(v)
    v_hat = [0 -v(3) v(2); 
             v(3) 0 -v(1);
             -v(2) v(1) 0];
end

%% Jacobian analytical
% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame, p, and
% magnet orientation with respect to world frame, R.
function J = J_analytical_sensor_B_to_world_pR(sens_pos, sens_or, magnet_pos, B_r, Volumn, R_star, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu_world_hat = R_star * [0;0;1];

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Calculate r in the world frame
    r_world = sens_pos - magnet_pos;
    r_world_hat = r_world/norm(r_world);
    
    % J of B with respect to p, evaluated at p*, R*
    J_position = -((3*B_r*Volumn)/(4*pi*norm(r_world)^4)) * ((eye(3)-5*(r_world_hat*r_world_hat.')) * ...
                 (r_world_hat.'*mu_world_hat) + ...
                 mu_world_hat*r_world_hat.' + r_world_hat*mu_world_hat.');

    % Second part of the jacobian, from sensor reading to magnet
    % orientation R
    
    % J of B with respect to R, evaluated at p*, R*
    e3 = [0;0;1];
    J_angle = (B_r*Volumn/(4*pi)) * (1/(norm(r_world)^3)) * (3*(r_world_hat*r_world_hat.')-eye(3)) * ...
              (-R_star*skew(e3));

    J = [J_position, J_angle];

    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end
