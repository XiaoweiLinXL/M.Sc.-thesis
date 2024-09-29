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

%% Genetic algorithm
lb = [-0.06/scale -0.06/scale ...
      -0.06/scale -0.06/scale ...
      -0.06/scale -0.06/scale ...
      -0.06/scale -0.06/scale ...
      4];
ub = [0.06/scale 0.06/scale ...
      0.06/scale 0.06/scale ...
      0.06/scale 0.06/scale ...
      0.06/scale 0.06/scale ...
      4];

initial_guess = [-0.06,-0.06,-0.06,0.06,0.06,-0.06,0.06,0.06, 4;
                 ];

options = optimoptions(@gamultiobj,'Display','iter', 'MaxStallGenerations', 500, ...
                       'MaxGenerations', 500, ...
                       'PopulationSize', 500, ...
                       'InitialPopulationMatrix', initial_guess);

fun = @(x) min_fun_orientation(x, magnet_conf, B_r, Volumn, type);
[sol, fval, exitflag, output] = gamultiobj(fun, length(lb), [], [], [], [], lb, ub, [], length(lb), options);

save('results_4_3_axis_multiobj_500gen_500pop_workspace012_largerworkspace')

%% Plot pareto front
load('results_9_3_axis_multiobj_1000gen_1000pop_workspace048_largerworkspace')
figure
scatter(-fval(:,1),-fval(:,2));
hold on
% scatter(0.012237, 1.62334e-7, '*')
grid on
xlabel('reciprocal condition number')
ylabel('min sigma_m')
% ylim([0.3e-7,1.8e-7])
title('Pareto front of the objectives')

%% Plot the configs on the border of the steps
obj_rcond = -fval(:,1);
obj_sigma_min = -fval(:,2);
index_1 = find(abs(obj_rcond-0.0405017) <= 1e-7);
index_2 = find(abs(obj_rcond-0.0551742) <= 1e-7);
index_3 = find(abs(obj_rcond-0.0565928) <= 1e-7);
index_4 = find(abs(obj_rcond-0.106318) <= 1e-6);
index_5 = find(abs(obj_rcond-0.107121) <= 1e-6);
index_6 = find(abs(obj_rcond-0.150436) <= 1e-6);
index_7 = find(abs(obj_rcond-0.151228) <= 1e-6);
index_8 = find(abs(obj_rcond-0.203725) <= 1e-6);

index_collection = [index_1,index_2,index_3,index_4,index_5,index_6,index_7,index_8];

% Store their corresponding sensor configs
sensor_config = [];
% Store the corresponding sensor configs
for i = 1:length(index_collection)
    sensor_config = [sensor_config; sol(index_collection(i),:)];
end

% Plot the sensor configs
for i = 1:size(sensor_config, 1)
    sensor_conf = sensor_config(i,:);
    sensor_num = sensor_conf(end);
    sensor_conf(end) = [];

    sensor_conf = reshape(sensor_conf, 2, []);
    sensor_conf = sensor_conf/10;
    
    % Example coordinates (replace with your own coordinates)
    x = sensor_conf(1,:);
    y = sensor_conf(2,:);
    
    % Calculate the radius using the average distance from the origin (0,0)
    distances = sqrt(x.^2 + y.^2);
    r = mean(distances);
    
    % Generate circle points for plotting
    theta = linspace(0, 2*pi, 100);
    circle_x = r * cos(theta);
    circle_y = r * sin(theta);
    
    % Plot the points and the circle
    figure;
    hold on;
    plot(x, y, 'bo', 'MarkerSize', 8);
    plot(circle_x, circle_y, '-','Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
    
    % Grid config
    scatter(-0.12, -0.12, 50, '*r')
    scatter(0.12, -0.12, 50, '*r')
    scatter(-0.12, 0.12, 50, '*r')
    scatter(0.12, 0.12, 50, '*r')
    scatter(0, 0.12, 50, '*r')
    scatter(0.12, 0, 50, '*r')
    scatter(0,-0.12,50,'*r')
    scatter(-0.12,0,50,'*r')
    scatter(0,0, 50, '*r')
    
    xlabel('X');
    ylabel('Y');
    xlim([-0.25 0.25])
    ylim([-0.25 0.25])
    axis equal;
    grid on;
    hold off;
    
    rcond_this_config = obj_rcond(index_collection(i));
    min_svd_this_config = obj_sigma_min(index_collection(i));
    title_str = sprintf('$\\chi = %f$, $\\sigma_{min} = %.1e$, Radius: %.2f', rcond_this_config, min_svd_this_config, r);
    title(title_str, 'Interpreter', 'latex')
end

%% Fit a circle to the configurations
sensor_conf = sensor_config(8,:);
sensor_num = sensor_conf(end);
sensor_conf(end) = [];

sensor_conf = reshape(sensor_conf, 2, []);
sensor_conf = sensor_conf/10;

% Example coordinates (replace with your own coordinates)
x = sensor_conf(1,:);
y = sensor_conf(2,:);

% Calculate the radius using the average distance from the origin (0,0)
distances = sqrt(x.^2 + y.^2);
r = mean(distances);

% Generate circle points for plotting
theta = linspace(0, 2*pi, 100);
circle_x = r * cos(theta);
circle_y = r * sin(theta);

% Plot the points and the circle
figure;
hold on;
plot(x, y, 'bo', 'MarkerSize', 8);
plot(circle_x, circle_y, '-','Color', [0, 0.4470, 0.7410], 'LineWidth', 2);

% Grid config
scatter(-0.12, -0.12, 50, '*r')
scatter(0.12, -0.12, 50, '*r')
scatter(-0.12, 0.12, 50, '*r')
scatter(0.12, 0.12, 50, '*r')
scatter(0, 0.12, 50, '*r')
scatter(0.12, 0, 50, '*r')
scatter(0,-0.12,50,'*r')
scatter(-0.12,0,50,'*r')
scatter(0,0, 50, '*r')

xlabel('X');
ylabel('Y');
xlim([-0.25 0.25])
ylim([-0.25 0.25])
% Specify the radius in the title
title(sprintf('Sensor location with fitted circle of Radius %.2f', r));
axis equal;
grid on;
hold off;



%% Chi distribution
number_of_sensor = 4;

axes = 3*number_of_sensor;
noise_std_one_axis = 5.3e-11;
noise_variance_one_axis = noise_std_one_axis^2;
noise_variance_measurement_pair = 2*noise_variance_one_axis;

x = 0:1e-12:0.3e-8;

% CDF of the norm's distribution
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

percentile998_noise_norm = threshold;

%%
resolution = 0.01;
figure(2)
scatter(-fval(:,1),-resolution*fval(:,2));
hold on
scatter(0.012237, resolution*1.62334e-7, '*')
grid on
xlabel('reciprocal condition number')
ylabel(["min signal's strength for resolution ", string(resolution)])
% ylim([0.3e-9,1.4e-9])
title('Pareto front of the objectives')

proportional_threshold = [];
for proportion = 1:4
    proportional_threshold = [proportional_threshold; proportion*threshold];
end
figure(2)
for i = 1:length(proportional_threshold)
    yline(proportional_threshold(i), '--');
end

legend('Pareto front', 'Grid result', ...
       '1 times 99.8 percentile of noise', ...
       '2 times 99.8 percentile of noise', ...
       '3 times 99.8 percentile of noise', ...
       '4 times 99.8 percentile of noise');

%% Extract the results that are better than the grid
load('results_4_3_axis_multiobj_1000gen_1000pop_workspace040_largerworkspace')
better_rcond_index = find(fval(:,1)<=-0.012237);
better_min_svd_index = find(fval(:,2)<=-1.62334e-7);
better_both_index = intersect(better_rcond_index, better_min_svd_index);

better_config = [];
% Store the corresponding sensor configs
for i = 1:length(better_both_index)
    better_config = [better_config; sol(better_both_index(i),:)];
end

% Plot the sensor configs
for i = 1:length(better_both_index)
    sensor_conf = better_config(i,:);
    sensor_num = sensor_conf(end);
    sensor_conf(end) = [];

    sensor_conf = reshape(sensor_conf, 2, []);
    sensor_conf = sensor_conf/10;
    
    % Better config
    figure
    for j = 1:size(sensor_conf,2)
        scatter(sensor_conf(1,j), sensor_conf(2,j), 100, 'blue')
        hold on
    end

    % Grid config
    scatter(-0.06, -0.06, 50, '*r')
    scatter(0.06, -0.06, 50, '*r')
    scatter(-0.06, 0.06, 50, '*r')
    scatter(0.06, 0.06, 50, '*r')

    xlim([-0.2,0.2])
    ylim([-0.2,0.2])
    grid on
    xlabel('x-coordinate')
    ylabel('y_coordinate')
    
    rcond = fval(better_both_index(i),1);
    min_svd = fval(better_both_index(i), 2);
%     title_str = sprintf('$\sigma_{min} = $%f, $\chi = $%f')
%     title('$\sigma_{min}$, $\chi$', 'Interpreter', 'latex')
end

%% Find 6 equal distance points on the pareto front based on sigma_min
load('results_9_3_axis_multiobj_1000gen_1000pop_workspace048_largerworkspace.mat')
rcond = -fval(:,1);
sigma_min = -fval(:,2);

% sorted_sigma_min = sort(sigma_min, 'descend');

% second_max_sigma_min = sorted_sigma_min(20);
max_sigma_min = max(sigma_min);
min_sigma_min = min(sigma_min);

num_points = 6;
interpolate_sigma_min = linspace(min_sigma_min, max_sigma_min, num_points); 

% Find the elements in the array that are closest to the interpolated
% values and their corresponding indices
closest_elements = zeros(1, num_points);
indices = zeros(1, num_points);
for i = 1:num_points
    [~, idx] = min(abs(sigma_min - interpolate_sigma_min(i)));
    closest_elements(i) = sigma_min(idx);
    indices(i) = idx;
end

% Store their corresponding sensor configs
sensor_config = [];
% Store the corresponding sensor configs
for i = 1:length(indices)
    sensor_config = [sensor_config; sol(indices(i),:)];
end

% Plot the sensor configs
for i = 1:length(indices)
    sensor_conf = sensor_config(i,:);
    sensor_num = sensor_conf(end);
    sensor_conf(end) = [];

    sensor_conf = reshape(sensor_conf, 2, []);
    sensor_conf = sensor_conf/10;
    
    % sensor config
    figure
    for j = 1:size(sensor_conf,2)
        scatter(sensor_conf(1,j), sensor_conf(2,j), 100, 'blue')
        hold on
    end

    % Grid config
    scatter(-0.12, -0.12, 50, '*r')
    scatter(0.12, -0.12, 50, '*r')
    scatter(-0.12, 0.12, 50, '*r')
    scatter(0.12, 0.12, 50, '*r')
    scatter(0, 0.12, 50, '*r')
    scatter(0.12, 0, 50, '*r')
    scatter(0,-0.12,50,'*r')
    scatter(-0.12,0,50,'*r')
    scatter(0,0, 50, '*r')

    xlim([-0.24,0.24])
    ylim([-0.24,0.24])
    grid on
    xlabel('x-coordinate')
    ylabel('y_coordinate')
    
    rcond_this_config = rcond(indices(i));
    min_svd_this_config = sigma_min(indices(i));
    title_str = sprintf('$\\chi = %f$, $\\sigma_{min} = %.1e$', rcond_this_config, min_svd_this_config);
    title(title_str, 'Interpreter', 'latex')
    axis equal
end

%% Range of required resolution affects pareto front
resolution = 0.01:0.01:1;
load('results_9_3_axis_multiobj_1000gen_1000pop_workspace048_largerworkspace')

% Find the max min sigma_min
max_min_sigma_min = -min(fval(:,2));

min_signal_strength_resolution = resolution * max_min_sigma_min;
ratio_to_noise = min_signal_strength_resolution/percentile998_noise_norm;

figure
plot(resolution, ratio_to_noise, 'LineWidth',2);
xlabel('resolution (dm)(rad)')
ylabel("ratio of signal's strength to noise");
grid on

%% Range of size of magnet affects pareto front
diameter = 3.175e-3:1e-3:5.08e-2;
load('results_4_3_axis_multiobj_1000gen_1000pop_workspace040_largerworkspace')

% Find the max min sigma_min
max_min_sigma_min = -min(fval(:,2));
volumn = (diameter/2).^3;
ratio_to_the_first_volumn = volumn/volumn(1);

min_signal_strength_size = max_min_sigma_min * ratio_to_the_first_volumn;

%%
signal_ratio_to_noise = [];
max_min_sigma_min = -min(fval(:,2)); % diameter 3.175e-3, resolution 1
for i = 1:length(resolution)
    for j=1:length(diameter)
        one_resolution = resolution(i);
        one_diameter = diameter(j);
        
        volumn_ratio_to_baseline = ((one_diameter/2)^3)/((diameter(1)/2)^3);

        min_signal_strength = volumn_ratio_to_baseline*one_resolution*max_min_sigma_min;
        ratio_to_noise = min_signal_strength/percentile998_noise_norm;
        signal_ratio_to_noise(i,j) = ratio_to_noise;
    end
end

[X, Y] = meshgrid(diameter, resolution);

figure
contour(X,Y,signal_ratio_to_noise);
colorbar;
xlabel("diameter")
ylabel("required resolution")
zlabel("")


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

%% Objective function
function obj = min_fun_orientation(sens_conf, magnet_conf, B_r, Volumn, type)
    
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    sens_conf = reshape(sens_conf, 2, []);
    
    % Add the line of zero as the z coordinate to the position of sensors
    z_coordinate = zeros(1, size(sens_conf, 2));
    sens_pos = [sens_conf; z_coordinate];

    % Default orientation
    default_or = [1;0;0;0];

    % Orientation for all sensors
    sens_or = repmat(default_or, 1, sens_num);
    
    % Collect the reciprocal condition number for each magnet configuration
    reciprocal_condition_number_set = [];
    min_svd_set = [];

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
        
        J = [];
        for i=1:sens_num
            J = [J;J_analytical_sensor_B_to_world_pR(sens_pos(:,i), sens_or(:,i), ...
                magnet_pos, B_r, Volumn, R_star, type)];
        end

        sigma = svd(J);

        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);

        reciprocal_condition_number_set = [reciprocal_condition_number_set; reciprocal_condition_number]; % Put the rcond for this magnet conf into the list
        min_svd_set = [min_svd_set; sigma(num_dof)];
    end
    
    % Minimize the negative of the min -> maximize the min
    obj = [-min(reciprocal_condition_number_set), -min(min_svd_set)];
end

function [c,ceq] = constraint(sens_conf, magnet_conf, B_r, Volumn, type)
    ceq = [];    % No nonlinear equalities at x.

    sens_num = 4;
    sens_conf = reshape(sens_conf, 2, []);
    
    % Add the line of zero as the z coordinate to the position of sensors
    z_coordinate = zeros(1, size(sens_conf, 2));
    sens_pos = [sens_conf; z_coordinate];

    % Default orientation
    default_or = [1;0;0;0];

    % Orientation for all sensors
    sens_or = repmat(default_or, 1, sens_num);
    
    % Collect the reciprocal condition number for each magnet configuration
    reciprocal_condition_number_set = [];
    min_svd_set = [];

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
        
        J = [];
        for i=1:sens_num
            J = [J;J_analytical_sensor_B_to_world_pR(sens_pos(:,i), sens_or(:,i), ...
                magnet_pos, B_r, Volumn, R_star, type)];
        end

        sigma = svd(J);

        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);

        reciprocal_condition_number_set = [reciprocal_condition_number_set; reciprocal_condition_number]; % Put the rcond for this magnet conf into the list
        min_svd_set = [min_svd_set; sigma(num_dof)];
    end

    min_svd = min(min_svd_set);

    c(1) = -min_svd + 1.6e-7;
end


