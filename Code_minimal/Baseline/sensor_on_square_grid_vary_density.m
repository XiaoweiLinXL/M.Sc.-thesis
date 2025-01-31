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
                                      % converted to corresponding unit
                                      % using the scale factor

B_r = 1.32; % Residual flux density (T)
Volumn = 24*(4/3)*pi*(norm(sph_dia)/2)^3; % volumn of the sphere dipole
                                       % in the unit of (specified unit)^3
                                       
sph_dip = B_r*Volumn/mu0; % spheres magnetic dipole mu_norm = B_r*V/mu0

mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole moment vector
mu_norm = norm(mu);

delta = 1e-7;
type = "3D";

%% Magnet configurations with orientation spherical
% phi_sphere = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4];
phi_sphere = linspace(0, 2*pi, 10);

theta_sphere = [pi/4, pi/2, 3*pi/4];
theta = [-pi/2, -pi/4, 0, pi/4, pi/2];
phi = [-pi/2, -pi/4, 0, pi/4, pi/2];
psi = [0];

[PHI_sphere, THETA_sphere, Theta, Phi, Psi] = ndgrid(phi_sphere, theta_sphere, theta, phi, psi);

PHI_sphere = reshape(PHI_sphere, [], 1);
THETA_sphere = reshape(THETA_sphere, [], 1);
Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);
Psi = reshape(Psi, [], 1);

sphere = [PHI_sphere, THETA_sphere, Theta, Phi, Psi].';

phi_sphere_extreme = [0];
theta_sphere_extreme = [0,pi];
[PHI_sphere, THETA_sphere, Theta, Phi, Psi] = ndgrid(phi_sphere_extreme, theta_sphere_extreme, theta, phi, psi);
PHI_sphere = reshape(PHI_sphere, [], 1);
THETA_sphere = reshape(THETA_sphere, [], 1);
Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);
Psi = reshape(Psi, [], 1);

sphere = [sphere, [PHI_sphere, THETA_sphere, Theta, Phi, Psi].'];

points = [];

r = 5*0.01; % 5*sqrt(3) cm
for i = 1:size(sphere, 2)
    phi = sphere(1, i);
    theta = sphere(2, i);
    x = r*sin(theta)*cos(phi)/scale;
    y = r*sin(theta)*sin(phi)/scale;
    z = (r*cos(theta)+0.15)/scale;

    points = [points, [x;y;z;sphere(3,i);sphere(4,i);sphere(5,i)]];
end

magnet_conf = points;


%% Sensor configuration
half_side_length = 0.16/scale;
LB = [-half_side_length, -half_side_length];
UB = [half_side_length, half_side_length];
z = 0;

number_of_sensor_on_side = 2:1:26;

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

% Remove the sensors in the middle
for i = 1:length(sensor_grid_config)
    sensor_config = sensor_grid_config{i};

    % Initialize a logical array for columns to keep
    keep_columns = true(1, size(sensor_config, 2));

    % Loop through each column to check the conditions
    for j = 1:size(sensor_config, 2)
        if abs(sensor_config(1, j)) ~= half_side_length && abs(sensor_config(2, j)) ~= half_side_length
            keep_columns(j) = false;  % Mark column for deletion
        end
    end

    % Keep only the columns that meet the condition
    sensor_grid_config{i} = sensor_config(:, keep_columns);
end

% two_sensor_conf = [-0.6,-0.6,0,1,0,0,0;-0.6,0.6,0,1,0,0,0].';
% three_sensor_conf = [-0.24,-0.24,0,1,0,0,0;0.24,-0.24,0,1,0,0,0;0,sqrt(2)*0.24,0,1,0,0,0].'/scale;

% five_sensor_conf = [-half_side_length,-half_side_length,0,1,0,0,0;...
%                     half_side_length,-half_side_length,0,1,0,0,0;...
%                     -half_side_length,half_side_length,0,1,0,0,0;...
%                     half_side_length,half_side_length,0,1,0,0,0;...
%                     0,0,0,1,0,0,0].';

% sensor_grid_config = [{three_sensor_conf}, sensor_grid_config];
% sensor_grid_config = [sensor_grid_config(1), {five_sensor_conf}, sensor_grid_config(2:end)];
% sensor_grid_config = [{five_sensor_conf}];
% sensor_grid_config= [{[0 0.2/scale 0 1 0 0 0; -0.2/scale -0.2/scale 0 1 0 0 0; 0.2/scale -0.2/scale 0 1 0 0 0].'}, sensor_grid_config];
% sensor_grid_config = [{[-0.2/scale 0 0 1 0 0 0; 0.2/scale 0 0 1 0 0 0].'},sensor_grid_config];

%% Sensors on circle
n = 4:100;

sensor_grid_config = {};
for i=1:length(n)
    theta = linspace(0, 2*pi, n(i)+1);
    theta(end) = [];

    x = half_side_length*sqrt(2) * cos(theta);
    y = half_side_length*sqrt(2) * sin(theta);
    conf = [x;y];
    conf = [conf;zeros(1,size(conf,2));ones(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2))];
    sensor_grid_config{i} = conf;
end

%%
% Loop over each sensor configuration
reciprocal_condition_number_all_sensor_config = [];
mean_reciprocal_condition_number = [];
min_singular_value_J_all_sensor_config = [];
% max_singular_value_J_all_sensor_config = [];
mean_singular_value_J_all_sensor_config = [];
mean_min_singular_value_J_all_sensor_config = [];
median_min_singular_value_J_all_sensor_config = [];

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
    max_singular_value_J = [];
    mean_singular_value_J = [];

    for j=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,j);
        theta = magnet_conf(4, j);
        phi = magnet_conf(5, j);
        psi = magnet_conf(6, j);

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

        J = [];
        for k = 1:sensor_number
            J = [J;J_analytical_sensor_B_to_world_pR(sensor_position(:,k), sensor_orientation(:,k), ...
                magnet_pos, B_r, Volumn, R_star, type)];
        end

        num_dof = 5;

        sigma_J = svd(J);

        rcond = sigma_J(num_dof)/sigma_J(1);
        min_svd_J = sigma_J(num_dof);
        max_svd_J = sigma_J(1);

        reciprocal_condition_number = [reciprocal_condition_number; rcond];
        min_singular_value_J = [min_singular_value_J; min_svd_J];
%         max_singular_value_J = [max_singular_value_J; max_svd_J];
        
        % Exclude the 0 sigma
        sigma_J(end) = [];
        mean_singular_value_J = [mean_singular_value_J; mean(sigma_J)];

    end

    % Put the minimum among all magnet configs into the sensor config
    % collection
    reciprocal_condition_number_all_sensor_config = [reciprocal_condition_number_all_sensor_config; min(reciprocal_condition_number)];
    mean_reciprocal_condition_number = [mean_reciprocal_condition_number; mean(reciprocal_condition_number)];
    min_singular_value_J_all_sensor_config = [min_singular_value_J_all_sensor_config; min(min_singular_value_J)];
%     max_singular_value_J_all_sensor_config = [max_singular_value_J_all_sensor_config; min(max_singular_value_J)];
    mean_singular_value_J_all_sensor_config = [mean_singular_value_J_all_sensor_config; mean(mean_singular_value_J)];
    mean_min_singular_value_J_all_sensor_config = [mean_min_singular_value_J_all_sensor_config; mean(min_singular_value_J)];
    median_min_singular_value_J_all_sensor_config = [median_min_singular_value_J_all_sensor_config; median(min_singular_value_J)];
end

%%
reciprocal_condition_number_all_sensor_config_3 = reciprocal_condition_number_all_sensor_config;
mean_reciprocal_condition_number_3 = mean_reciprocal_condition_number;
min_singular_value_J_all_sensor_config_3 = min_singular_value_J_all_sensor_config;
mean_singular_value_J_all_sensor_config_3 = mean_singular_value_J_all_sensor_config;
mean_min_singular_value_J_all_sensor_config_3 = mean_min_singular_value_J_all_sensor_config;
median_min_singular_value_J_all_sensor_config_3 = median_min_singular_value_J_all_sensor_config;

%%
number_of_sensor = [[2:1:10].^2];
% number_of_sensor = [3,number_of_sensor];
% number_of_sensor = [number_of_sensor(1),5,number_of_sensor(2:end)];
number_of_sensor_peri = 4:4:100;
number_of_sensor_circle = 4:1:100;

figure(1)
plot(number_of_sensor, reciprocal_condition_number_all_sensor_config_1, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2)
hold on
% plot(number_of_sensor_peri, reciprocal_condition_number_all_sensor_config_2, '-', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2)
% hold on 
% plot(number_of_sensor_circle, reciprocal_condition_number_all_sensor_config_3, '-', 'Color', [0.9290, 0.6940, 0.1250], "linewidth", 2)
% hold on 
% plot(number_of_sensor_peri, reciprocal_condition_number_all_sensor_config_4, '--', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2)
% hold on 
% plot(number_of_sensor, reciprocal_condition_number_all_sensor_config_5, '-', 'Color', [0.9290, 0.6940, 0.1250], "linewidth", 2)
% hold on 
% plot(number_of_sensor, reciprocal_condition_number_all_sensor_config_6, '--', 'Color', [0.9290, 0.6940, 0.1250], "linewidth", 2)
% figure(3)
% plot(number_of_sensor, mean_reciprocal_condition_number_1)

xlabel('Number of sensors' , 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylim([0,0.3])
% title('Min reciprocal condition number vs number of sensor')
grid on
% legend('grid 0.28m', 'square 0.24m, adding only on the sides', ...
%        'grid 0.48m', 'square 0.48m, adding only on the sides')

legend('grid 0.32m', 'circle 0.5m', 'circle 0.32m')

% c = 0.01; % Required accuracy
figure(2)
plot(number_of_sensor, mean_min_singular_value_J_all_sensor_config_1, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2)
hold on
plot(number_of_sensor_peri, mean_min_singular_value_J_all_sensor_config_2, '-', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2)
hold on
plot(number_of_sensor_circle, median_min_singular_value_J_all_sensor_config_3, '-', 'Color', [0.9290, 0.6940, 0.1250], "linewidth", 2)
hold on
% plot(number_of_sensor_peri, 0.01*min_singular_value_J_all_sensor_config_4, '--', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2)
% hold on
% plot(number_of_sensor, 0.01*min_singular_value_J_all_sensor_config_5, '-', 'Color', [0.9290, 0.6940, 0.1250], "linewidth", 2)
% hold on
% plot(number_of_sensor, 0.01*min_singular_value_J_all_sensor_config_6, '--', 'Color', [0.9290, 0.6940, 0.1250], "linewidth",2)
legend('grid 0.32m', 'square 0.32m', 'circle 0.32m')

xlabel('Number of sensors', 'FontSize', 12 , 'FontWeight', 'bold' ,'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
% ylim([0,9e-6])
% title('Min sigma_m of J vs number of sensor')
grid on

% figure(3)
% plot(number_of_sensor, mean_min_singular_value_J_all_sensor_config_1, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2)

%%
number_of_sensor = [[2:1:10].^2];

figure
plot(number_of_sensor, reciprocal_condition_number_all_sensor_config_3)
figure
plot(number_of_sensor, min_singular_value_J_all_sensor_config_3)
figure
plot(number_of_sensor, mean_singular_value_J_all_sensor_config_3)

%% Chi distribution
number_of_sensor_range = 3:1:100;
percentile998_noise_norm = [];

for index = 1:length(number_of_sensor_range)
    number_of_sensor = number_of_sensor_range(index);
    axes = 3*number_of_sensor;
    noise_std_one_axis = 5e-08;
    noise_variance_one_axis = noise_std_one_axis^2;
    noise_variance_measurement_pair = 2*noise_variance_one_axis;
    
    x = 0:1e-10:1e-6;
    
    cdf_analytical = [];
    for i = 1:length(x)
        cdf_analytical_onepoint = gammainc(axes/2,(x(i)^2)/(2*noise_variance_measurement_pair));
        cdf_analytical = [cdf_analytical, cdf_analytical_onepoint];
    end
    cdf_analytical=1-cdf_analytical;
    
    % figure
    % plot(x, cdf_analytical);
        
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
number_of_sensor = [number_of_sensor_on_side.^2];
% number_of_sensor = [2,3,number_of_sensor];
% plot(number_of_sensor, min_singular_value_J_all_sensor_config*required_accuracy, "LineWidth", 2)
% hold on
plot(number_of_sensor_range, percentile998_noise_norm*ratio, '-', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth',2)
legend('grid 0.34m', 'square 0.34m', 'circle 0.34m', '99.8 percentile noise')


% grid on
% xlabel("number of sensor")
% ylabel("99.8 percentile of the norm of the noise vector / min change in B")
% title("Min change in B vs Noise vector's norm")
% legend("Minimum change in B, 1mm-0.5deg, accuracy", strcat(string(ratio), " times 99.8 percentile of the noise norm"))

%% Find out the intersection points
min_singular_value_J_all_sensor_config_interpolated = interp1(number_of_sensor, ...
    min_singular_value_J_all_sensor_config, number_of_sensor_range, "spline");

figure(4)
plot(number_of_sensor_range, min_singular_value_J_all_sensor_config_interpolated*required_accuracy, "LineWidth", 2)
hold on
plot(number_of_sensor_range, percentile998_noise_norm*ratio, "LineWidth", 2)
grid on
xlabel("number of sensor")
ylabel("99.8 percentile of the norm of the noise vector / min change in B")
title("Min change in B vs Noise vector's norm")
legend("Minimum change in B, 1mm-0.5deg, accuracy", strcat(string(ratio), " times 99.8 percentile of the noise norm"))


for accuracy = 0.00227:0.000001:0.00284
    min_B_for_one_accuracy = min_singular_value_J_all_sensor_config_interpolated*accuracy;
%     figure
%     plot(number_of_sensor_range, min_B_for_one_accuracy, "LineWidth", 2)
%     hold on
%     plot(number_of_sensor_range, percentile998_noise_norm*ratio, "LineWidth", 2)

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


