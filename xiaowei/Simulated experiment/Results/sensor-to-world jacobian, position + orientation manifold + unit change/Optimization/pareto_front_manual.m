%%
close all
clear
clc

%% Spheres
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
phi_sphere = linspace(0, 2*pi, 9);

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

%% Experiment configuration
theta = linspace(0, 2*pi, 9);
theta(end) = [];
phi = linspace(0, pi, 5);
phi(end) = [];
psi = [0];

[Theta, Phi, Psi] = ndgrid(theta, phi, psi);

Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);
Psi = reshape(Psi, [], 1);

Orientation = [Theta, Phi, Psi].';

Position = [0, 0, 0.1; 0, 0, 0.15; 0, 0, 0.2; 0.05, 0, 0.15; -0.05, 0, 0.15; 0, 0.05, 0.15; 0, -0.05, 0.15].'/scale;
Position = [0.0;0.0;0.15]/scale;

nPos = size(Position, 2); % Number of positions (7)
nOri = size(Orientation, 2); % Number of orientations (64)

% Initialize the new array
magnet_conf = zeros(6, nPos * nOri);

% Fill the new array
index = 1;
for i = 1:nPos
    for j = 1:nOri
        % Combine position i and orientation j
        magnet_conf(:, index) = [Position(:, i); Orientation(:, j)];
        index = index + 1;
    end
end

magnet_conf = magnet_conf(:,1);


%% Define the 9 sensors configuration located on the circle of different radius
% Number of sensor
sens_num = 9;

% Define the radius
radius = 0.01:0.001:1;

% Initialize a cell array to hold all coordinates for each value of k
all_coordinates = cell(1, length(radius));

for i = 1:length(radius)
    % Angles for the points on the circle
    angles = linspace(0, 2*pi, sens_num+1);
    angles(end) = []; % Remove the last element to get k points
    
    % Calculate coordinates for the circle
    x = radius(i) * cos(angles);
    y = radius(i) * sin(angles);
    
    % Store the coordinates
    all_coordinates{i} = [x', y'].';
end

% Iterate throught the sensor configs to get the pareto front
rcond_all_sens_conf_circle = [];
average_rcond_all_sens_conf_circle = [];
median_rcond_all_sens_conf_circle = [];
min_svd_all_sens_conf_circle = [];
average_svd_all_sens_conf_circle = [];
median_min_svd_all_sens_conf_circle = [];
mean_min_svd_all_sens_conf_circle = [];

for i = 1:length(all_coordinates)
    sens_pos = all_coordinates{i};
    sens_num = size(sens_pos,2);

%     % Add a sensor in the middle
%     sens_pos = [sens_pos, [0;0]];
%     sens_num = sens_num+1;

    % Add the line of zero as the z coordinate to the position of sensors
    z_coordinate = zeros(1, size(sens_pos, 2));
    sens_pos = [sens_pos; z_coordinate];
    sens_pos = sens_pos/scale;

    % sens_pos = [0.0463,-0.0012,0.0015;
    %             -0.0039,0.0507,0.0038;
    %             -0.0537,-0.0015,0.0025;
    %             -0.0040,-0.0515,0.0011].'/scale;

    % sens_pos = [0.0972,-0.0013,0.0002;
    %             -0.0039,0.1015,0.0030;
    %             -0.1054,-0.0017,0.0029;
    %             -0.0038,-0.1017,-0.0003].'/scale;



    % Default orientation
    default_or = [1;0;0;0];

    % Orientation for all sensors
    sens_or_unitary = repmat(default_or, 1, sens_num);
    
    rcond_one_sens_conf = [];
    min_svd_one_sens_conf = [];
    average_sigma_one_sens_conf = [];

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
    
        % R_star = Rx_1 * Ry * Rx_2;
        R_star = Ry * Rx_1;
    
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

        % Exclude the 0 sigma_min
        sigma_scaled(end) = [];
        average_sigma_one_sens_conf = [average_sigma_one_sens_conf; mean(sigma_scaled)];
    end

    index_min_rcond = find(rcond_one_sens_conf == min(rcond_one_sens_conf));
    index_min_minsvd = find(min_svd_one_sens_conf == min(min_svd_one_sens_conf));

    % Put the min value in the workspace into the list
    rcond_all_sens_conf_circle = [rcond_all_sens_conf_circle; min(rcond_one_sens_conf)];
    average_rcond_all_sens_conf_circle = [average_rcond_all_sens_conf_circle; mean(rcond_one_sens_conf)];
    median_rcond_all_sens_conf_circle = [median_rcond_all_sens_conf_circle; median(rcond_one_sens_conf)];
    min_svd_all_sens_conf_circle = [min_svd_all_sens_conf_circle; min(min_svd_one_sens_conf)];
    mean_min_svd_all_sens_conf_circle = [mean_min_svd_all_sens_conf_circle; mean(min_svd_one_sens_conf)];
    median_min_svd_all_sens_conf_circle = [median_min_svd_all_sens_conf_circle; median(min_svd_one_sens_conf)];
    average_svd_all_sens_conf_circle = [average_svd_all_sens_conf_circle; mean(average_sigma_one_sens_conf)];
end

%% Find the pareto front
% figure
% scatter(rcond_all_sens_conf_circle, min_svd_all_sens_conf_circle)
% figure
% scatter(rcond_all_sens_conf_circle, average_svd_all_sens_conf_circle)
% hold on
% % load('results_16_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma.mat')
% % scatter(-fval(:,1),-fval(:,2)*24, 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
% 
% figure
% scatter(average_rcond_all_sens_conf_circle, average_svd_all_sens_conf_circle)
% % hold on
% % load('results_4_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_mean_rcond.mat')
% % scatter(-fval(:,1),-fval(:,2)*24, 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
% 
figure
scatter(rcond_all_sens_conf_circle, median_min_svd_all_sens_conf_circle)
hold on
load("results_5_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_median_sigma_min.mat")
scatter(-fval(:,1),-fval(:,2), 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
xlabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$median_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')

% % figure
% scatter(median_rcond_all_sens_conf_circle, median_min_svd_all_sens_conf_circle)

figure
hold on
plot(radius, min_svd_all_sens_conf_circle, LineWidth=2)
plot(radius, mean_min_svd_all_sens_conf_circle, LineWidth=2)
plot(radius, average_svd_all_sens_conf_circle, LineWidth=2)
plot(radius, median_min_svd_all_sens_conf_circle, LineWidth=2)
legend({'$min \ \sigma_{min}$',  '$mean \ \sigma_{min}$', '$mean \ \sigma_{mean}$'}, 'Interpreter', 'latex')

%% Define the 9 sensors configuration located on the grid of different size
side_length = 0.01:0.001:1;

rcond_all_sens_conf_grid = [];
min_svd_all_sens_conf_grid = [];
for i = 1:length(side_length)
    
    rows = 2;
    columns = 2;
    sens_num = rows*columns;
    coordinate_cell = cell(rows,columns);
    
    for row = 1:rows
        for column = 1:columns
            coordinate_cell{row, column} = [row-1, column-1]; % Store the coordinates as a vector (i-1, j-1)
        end
    end
    
    range = 2*side_length(i);
    
    for row = 1:rows
        for column = 1:columns
            coords = coordinate_cell{row,column};
            orig_x = coords(1);
            orig_y = coords(2);
    
            scaled_x = orig_x * (range / (rows-1)) - side_length(i);
            scaled_y = orig_y * (range / (columns-1)) - side_length(i);
    
            coordinate_cell{row,column} = [scaled_x, scaled_y];
        end
    end
    
    % Initialize the 2x(rows*columns) array
    sens_pos = zeros(2, rows * columns);
    
    % Fill the 2x(rows*columns) array with the scaled coordinates
    index = 1;
    for row = 1:rows
        for column = 1:columns
            sens_pos(:, index) = coordinate_cell{row, column};
            index = index + 1;
        end
    end

    % % Add a sensor in the middle
    % sens_pos = [sens_pos, [0;0]];
    % sens_num = sens_num+1;

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
    rcond_all_sens_conf_grid = [rcond_all_sens_conf_grid; min(rcond_one_sens_conf)];
    min_svd_all_sens_conf_grid = [min_svd_all_sens_conf_grid; min(min_svd_one_sens_conf)];
end

%%
% figure
scatter(rcond_all_sens_conf_grid, min_svd_all_sens_conf_grid)
hold on

%%
% Define the side length of the square
L = 0.01:0.001:1; % Example side length

rcond_all_sens_conf_square = [];
min_svd_all_sens_conf_square = [];

for i = 1:length(L)
    side_length = L(i);
    % Calculate the distance between points
    sens_num = 4;
    distance_between_points = (4 * 2 * side_length) / sens_num;
    
    % Initialize arrays to hold the coordinates
    x_coords = zeros(sens_num, 1);
    y_coords = zeros(sens_num, 1);
    
    % Place the first point at the bottom-left corner of the centered square
    x_coords(1) = -side_length;
    y_coords(1) = -side_length;
    
    % Calculate the positions of the remaining points
    for i = 2:sens_num
        % Calculate the position along the perimeter
        pos = (i-1) * distance_between_points;
        
        if pos <= 2 * side_length
            % Bottom edge
            x_coords(i) = -side_length + pos;
            y_coords(i) = -side_length;
        elseif pos <= 4 * side_length
            % Right edge
            x_coords(i) = side_length;
            y_coords(i) = -side_length + (pos - 2 * side_length);
        elseif pos <= 6 * side_length
            % Top edge
            x_coords(i) = side_length - (pos - 4 * side_length);
            y_coords(i) = side_length;
        else
            % Left edge
            x_coords(i) = -side_length;
            y_coords(i) = side_length - (pos - 6 * side_length);
        end
    end
    
    sens_pos = [x_coords.'; y_coords.'];

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
    rcond_all_sens_conf_square = [rcond_all_sens_conf_square; min(rcond_one_sens_conf)];
    min_svd_all_sens_conf_square = [min_svd_all_sens_conf_square; min(min_svd_one_sens_conf)];
end

%%
% figure
scatter(rcond_all_sens_conf_square, min_svd_all_sens_conf_square)
%% Pareto front
figure
hold on
plot_points_with_arrows(rcond_all_sens_conf_circle, median_min_svd_all_sens_conf_circle, [0.8500, 0.3250, 0.0980])
hold on
% plot_points_with_arrows(rcond_all_sens_conf_grid, min_svd_all_sens_conf_grid, [0.4940, 0.1840, 0.5560])
% hold on
% plot_points_with_arrows(rcond_all_sens_conf_square, min_svd_all_sens_conf_square, [0.6350, 0.0780, 0.1840])
% hold on
 
grid on

load('results_9_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_median_sigma_min.mat')
scatter(-fval(:,1),-fval(:,2)*24, 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
hold on
xlabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
xlim([0 0.30])
ylim([0 1.8e-5])

legend('sensors on circle', 'pareto front form GA', 'Location', 'southeast')

%% Chi distribution
number_of_sensor_range = 4:1:4;
percentile998_noise_norm = [];

for index = 1:length(number_of_sensor_range)
    number_of_sensor = number_of_sensor_range(index);
    axes = 3*number_of_sensor;
    noise_std_one_axis = 5e-08;
    noise_variance_one_axis = noise_std_one_axis^2;
    noise_variance_measurement_pair = 2*noise_variance_one_axis;
    
    x = 0:1e-10:5e-6;
    
    cdf_analytical = [];
    for i = 1:length(x)
        cdf_analytical_onepoint = gammainc(axes/2,(x(i)^2)/(2*noise_variance_measurement_pair));
        cdf_analytical = [cdf_analytical, cdf_analytical_onepoint];
    end
    cdf_analytical=1-cdf_analytical;
        
    for j = 1:length(x)
        if cdf_analytical(j) >= 0.99
            threshold = x(j);
            break
        end
    end

    percentile998_noise_norm = [percentile998_noise_norm, threshold];
end

% Cylinder volumn
face_dia = 6.35e-3*[0; 0; 1] / scale; % face diameter 2mm
                                      % converted to corresponding unit
                                      % using the scale factor
magnet_length = 6.35*2e-3 / scale;                                  

Volumn_cylinder = pi*(norm(face_dia)/2)^2*magnet_length; % volumn of the sphere dipole
                                       % in the unit of (specified unit)^3

% Ratio between the volumn of the cylinder and the sphere
volumn_ratio = Volumn_cylinder/Volumn;

% Scale the noise level using to the required accuracy and volumn
accuracy = 0.5;
percentile998_noise_norm_accuracy = percentile998_noise_norm / (accuracy*volumn_ratio);
percentile998_noise_norm_accuracy_collection = [1,2,3,30]*percentile998_noise_norm_accuracy;

for i = 1:length(percentile998_noise_norm_accuracy_collection)
    hline = plot(xlim, [percentile998_noise_norm_accuracy_collection(i) percentile998_noise_norm_accuracy_collection(i)], '-', 'LineWidth', 1.5, 'HandleVisibility','off');
    hline.Color = [0, 0.4470, 0.7410]; % [0, 0.4470, 0.7410]
    hline.Color(4) = 0.2;
end

yyaxis right
h = gca;
h.YColor = [0, 0.4470, 0.7410];
ylim([0 1.0e-7])
yticks(percentile998_noise_norm_accuracy_collection); % Define the positions of the ticks
yticklabels({'1x', '2x', '3x', '30x'}); 
ylabel('SNR ($d=h=6.35$mm)', 'Interpreter', 'latex', 'Color', 'black', 'FontSize',11)

accuracy = 0.01:0.01:0.1; % 1mm-10mm, 0.5deg-5.0deg
for i = 1:length(accuracy)
    accuracy_one = accuracy(i);
    
end


%% Plateau plot of gain margin
load('results_9_3_axis_multiobj_1000gen_1000pop_workspace048_largerworkspace.mat')
max_min_svd = max(-fval(:,2));
min_min_svd = min(-fval(:,2));
accuracy = 0.001:0.001:0.1; % 0.1mm to 1cm
volumn = pi*(norm(2e-3/scale)/2)^2*(2e-3/scale):1e-6:pi*(norm(5e-3/scale)/2)^2*(5e-3/scale); % 2mm to 5mm cylinder

gain_margin = zeros(length(accuracy), length(volumn));
gain_margin_min = zeros(length(accuracy), length(volumn));

for i = 1:length(accuracy) 
    for j = 1:length(volumn)
        volumn_ratio = volumn(j)/((4/3)*pi*(norm(3.175e-3/scale)/2)^3);
        percentile998_noise_norm_accuracy = percentile998_noise_norm / (accuracy(i)*volumn_ratio);
        max_snr = max_min_svd / percentile998_noise_norm_accuracy;
        min_snr = min_min_svd / percentile998_noise_norm_accuracy;
        gain_margin(i,j) = max_snr;
        gain_margin_min(i,j) = min_snr;
    end
end

figure
imagesc(volumn, accuracy, gain_margin);
set(gca, 'YDir', 'normal'); % To set the y-axis direction to normal (instead of reverse)
colorbar; % Add a color bar to indicate the color scale

figure
imagesc(volumn, accuracy, gain_margin_min);
set(gca, 'YDir', 'normal'); % To set the y-axis direction to normal (instead of reverse)
colorbar; % Add a color bar to indicate the color scale

rcond_sorted = sort(-fval(:,1));
min_svd_sorted = sort(-fval(:,2), 'descend');

rcond_margin = zeros(length(accuracy), length(volumn));
for i = 1:size(gain_margin_min,1)
    for j = 1:size(gain_margin_min,2)
        if gain_margin_min(i, j) < 1
            volumn_ratio = volumn(j)/((4/3)*pi*(norm(3.175e-3/scale)/2)^3);
            percentile998_noise_norm_accuracy = percentile998_noise_norm / (accuracy(i)*volumn_ratio);
            index = find(min_svd_sorted > percentile998_noise_norm_accuracy);

            if ~isempty(index)
                rcond_margin(i, j) = rcond_sorted(max(index));
            end
        else
            rcond_margin(i, j) = max(rcond_sorted);
        end
    end
end

figure
imagesc(volumn, accuracy, rcond_margin);
set(gca, 'YDir', 'normal'); % To set the y-axis direction to normal (instead of reverse)
colorbar; % Add a color bar to indicate the color scale

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

%% Plot points with arrows
function [] = plot_points_with_arrows(x_coord, y_coord, color)
    side_length = 2 * 0.01:0.001:1;
    % Define the points
    points = [x_coord, y_coord];
    
    % Plot the points (for reference)
    plot(points(:,1), points(:,2), '-', 'Color',color, 'MarkerSize',2);
    hold on;
    
    % Set axis limits for better visualization
    xlim([0 0.30]);
    ylim([0 1.8e-5]);
    
    % Get the current axes position
    ax = gca;
    ax.Units = 'normalized';
    ax_pos = ax.Position;
    
    % Normalize the coordinates to the figure
    x_range = ax.XLim(2) - ax.XLim(1);
    y_range = ax.YLim(2) - ax.YLim(1);
    
    % Normalize each point
    norm_points = zeros(size(points));
    for i = 1:size(points, 1)
        norm_points(i, 1) = (points(i, 1) - ax.XLim(1)) / x_range * ax_pos(3) + ax_pos(1);
        norm_points(i, 2) = (points(i, 2) - ax.YLim(1)) / y_range * ax_pos(4) + ax_pos(2);
    end
    
    % Add arrows at equal distances
    num_arrows = 10;
    step = 1;
    total_points = size(points, 1);
    interval = floor((total_points - 20) / (num_arrows + 1)); % Adjust interval calculation

    for i = 0:num_arrows-1
        start_idx = 20 + i * interval; % Start from the 10th coordinate
        end_idx = start_idx + step;
        if end_idx <= total_points
            arrow = annotation('arrow', [norm_points(start_idx, 1), norm_points(end_idx, 1)], ...
                                        [norm_points(start_idx, 2), norm_points(end_idx, 2)], ...
                                        'Color', color, 'LineWidth', 2); % Set the color and width of the arrow
        end
    end

    
    % Labels
    xlabel('X-axis');
    ylabel('Y-axis');
    hold off;
end