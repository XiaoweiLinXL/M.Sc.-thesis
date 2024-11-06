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


%% For a certain inner/outer radius, change the density of the sensor grid
% Define the base radius and maximum radius
base_radius = 0.028; % You can change this to any value for the base radius
max_radius = 5.5 * base_radius; % You can change the multiplier to set the max radius

circular_grid_config = {};
for i = 2:10
    % Number of points
    n = i;
    
    % Number of circles
    num_circles = i;
    
    % Calculate the radii for the circles
    radii = linspace(base_radius, max_radius, num_circles);
    
    % Initialize arrays to hold all coordinates
    all_coordinates = cell(1, num_circles);
    
    figure;
    hold on;
    
    % Loop through each circle
    for k = 1:num_circles
        % Current radius
        radius = radii(k);
        
        % Angles for the points on the circle
        if k == 1
            angles = linspace(0, 2*pi, n+1);
        else
            angles = linspace(0, 2*pi, n+1) + (pi/n) * (k - 1);
        end
        angles(end) = []; % Remove the last element to get n points
        
        % Calculate coordinates for the circle
        x = radius * cos(angles);
        y = radius * sin(angles);
        
        % Store the coordinates
        all_coordinates{k} = [x', y'];
        
        % Plot the circle and points
        theta = linspace(0, 2*pi, 100);
        circle_x = radius * cos(theta);
        circle_y = radius * sin(theta);
        
        plot(circle_x, circle_y, 'b'); % Plot the circle
        plot(x, y, 'ro'); % Plot the points on the circle
    end
    
    axis equal;
    title('Equally Distant Points on Multiple Circles');
    xlabel('x');
    ylabel('y');
    grid on;
    
    % Concatenate the coordinates for all circles
    concatenated_coordinates = [];
    for i = 1:num_circles
        concatenated_coordinates = [concatenated_coordinates; all_coordinates{i}];
    end
    
    % Display the coordinates for all circles
    disp('Concatenated coordinates for all circles:');
    disp(concatenated_coordinates);

    circular_grid_config{end+1} = concatenated_coordinates.';
end

three_sensor_pos = [max_radius,0;-max_radius/2, max_radius*sqrt(3)/2 ;  -max_radius/2,-max_radius*sqrt(3)/2].';

circular_grid_config = [{three_sensor_pos}, circular_grid_config];

%% Iterate through the sensor config collection
rcond_all_sens_conf = [];
min_svd_all_sens_conf = [];

for i = 1:length(circular_grid_config)
    sensor_config_one = circular_grid_config{i};
    sens_num = size(sensor_config_one,2);
    sens_pos = sensor_config_one;

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
    end
    
    % Put the min value in the workspace into the list
    rcond_all_sens_conf = [rcond_all_sens_conf; min(rcond_one_sens_conf)];
    min_svd_all_sens_conf = [min_svd_all_sens_conf; min(min_svd_one_sens_conf)];
end

%%
rcond_all_sens_conf_2 = rcond_all_sens_conf;
min_svd_all_sens_conf_2 = min_svd_all_sens_conf;

%% Plots
sensor_number = [3,4,9,16,25,36,49,64,81,100];
figure(1)
plot(sensor_number, rcond_all_sens_conf_1, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2);
hold on
plot(sensor_number, rcond_all_sens_conf_2, 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2);
hold on
% ylim([0 1])
xlim([0 100])
xlabel("number of sensor")
ylabel("Minimum reciprocal condition number in the workspace")
title("Minimum reciprocal condition number in the workspace vs number of sensor with fixed radius")
grid on

figure(2)
plot(sensor_number, 0.01*min_svd_all_sens_conf_1, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2);
hold on
plot(sensor_number, 0.01*min_svd_all_sens_conf_2, '-', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2);
hold on
% ylim([0 2.5e-7])
xlim([0 100])
xlabel("number of sensor")
ylabel("Minimum of the \sigma_{min} in the workspace")
title("Minimum of the \sigma_{min} in the workspace vs number of sensor with fixed radius")
grid on



%% Adding sensor only on the edge
% Define the radius
radius = 0.154;

% Initialize a cell array to hold all coordinates for each value of k
all_coordinates = cell(98, 1);

% Loop through the number of points from 3 to 100
for k = 3:100
    % Angles for the points on the circle
    angles = linspace(0, 2*pi, k+1);
    angles(end) = []; % Remove the last element to get k points
    
    % Calculate coordinates for the circle
    x = radius * cos(angles);
    y = radius * sin(angles);
    
    % Store the coordinates
    all_coordinates{k-2} = [x', y'].';
end

%% Iterate through the sensor config collection
rcond_all_sens_conf = [];
min_svd_all_sens_conf = [];

for i = 1:length(all_coordinates)
    sensor_config_one = all_coordinates{i};
    sens_num = size(sensor_config_one,2);
    sens_pos = sensor_config_one;

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
    end
    
    % Put the min value in the workspace into the list
    rcond_all_sens_conf = [rcond_all_sens_conf; min(rcond_one_sens_conf)];
    min_svd_all_sens_conf = [min_svd_all_sens_conf; min(min_svd_one_sens_conf)];
end

%%
rcond_all_sens_conf_4 = rcond_all_sens_conf;
min_svd_all_sens_conf_4 = min_svd_all_sens_conf;

%% Plots
sensor_number = 3:1:100;
figure(1)
plot(sensor_number, rcond_all_sens_conf_3, '--', 'Color', [0, 0.4470, 0.7410], "linewidth", 2);
hold on
plot(sensor_number, rcond_all_sens_conf_4, '--', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2);
ylim([0 0.3])
xlim([0 100])
xlabel("number of sensor")
ylabel("Minimum reciprocal condition number in the workspace")
title("Minimum reciprocal condition number in the workspace vs number of sensor with fixed radius")
grid on
legend('circular grid with radius 0.165m', 'circular grid with radius 0.154m', 'circular grid with radius 0.165m, adding only on the periphery',...
    'circular grid with radius 0.154m, adding only on the periphery')

figure(2)
plot(sensor_number, 0.01*min_svd_all_sens_conf_3, '--', 'Color', [0, 0.4470, 0.7410], "linewidth", 2);
hold on
plot(sensor_number, 0.01*min_svd_all_sens_conf_4, '--', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2);
hold on
% ylim([0 2.5e-7])
xlim([0 100])
xlabel("number of sensor")
ylabel("Minimum of the \sigma_{min} in the workspace")
title("Minimum of the \sigma_{min} in the workspace vs number of sensor with fixed radius")
grid on


%% Chi distribution
number_of_sensor_range = 3:1:100;
percentile998_noise_norm = [];

for index = 1:length(number_of_sensor_range)
    number_of_sensor = number_of_sensor_range(index);
    axes = 3*number_of_sensor;
    noise_std_one_axis = 5.3e-11;
    noise_variance_one_axis = noise_std_one_axis^2;
    noise_variance_measurement_pair = 2*noise_variance_one_axis;
    
    x = 0:1e-12:0.3e-8;
    
    cdf_analytical = [];
    for i = 1:length(x)
        cdf_analytical_onepoint = gammainc(axes/2,(x(i)^2)/(2*noise_variance_measurement_pair));
        cdf_analytical = [cdf_analytical, cdf_analytical_onepoint];
    end
    cdf_analytical=1-cdf_analytical;
        
    for j = 1:length(x)
        if cdf_analytical(j) >= 0.998
            threshold = x(j);
            break
        end
    end

    percentile998_noise_norm = [percentile998_noise_norm, threshold];

end

ratio = 1;
required_accuracy = 0.01;
figure(2)
% number_of_sensor = [2,3,number_of_sensor];
% plot(number_of_sensor, min_singular_value_J_all_sensor_config*required_accuracy, "LineWidth", 2)
% hold on
plot(number_of_sensor_range, percentile998_noise_norm*ratio, '-', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth',2)
legend('circular grid with radius 0.165m', ...
       'circular grid with radius 0.154m', ...
       'circular grid with radius 0.165m, adding only on the periphery', ....
       'circular grid with radius 0.154m, adding only on the periphery', ...
       'noise')

% grid on
% xlabel("number of sensor")
% ylabel("99.8 percentile of the norm of the noise vector / min change in B")
% title("Min change in B vs Noise vector's norm")
% legend("Minimum change in B, 1mm-0.5deg, accuracy", strcat(string(ratio), " times 99.8 percentile of the noise norm"))

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
