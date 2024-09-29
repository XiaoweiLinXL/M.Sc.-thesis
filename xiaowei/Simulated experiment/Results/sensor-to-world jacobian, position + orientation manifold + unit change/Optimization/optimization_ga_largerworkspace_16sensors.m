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

%% Genetic algorithm
lb = [-0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      -0.5/scale -0.5/scale ...
      16];
ub = [0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      0.5/scale 0.5/scale ...
      16];

angles = linspace(0, 2*pi, 16+1);
angles(end) = [];
x = 0.12 * cos(angles);
y = 0.12 * sin(angles);
initial_guess = [x; y];
initial_guess = reshape(initial_guess, [], 1);
initial_guess = [initial_guess.', 16];

options = optimoptions(@gamultiobj,'Display','iter', 'MaxStallGenerations', 1000, ...
                       'MaxGenerations', 1000, ...
                       'PopulationSize', 1000, ...
                       'InitialPopulationMatrix', initial_guess);

fun = @(x) min_fun_orientation(x, magnet_conf, B_r, Volumn, type);
[sol, fval, exitflag, output] = gamultiobj(fun, length(lb), [], [], [], [], lb, ub, [], length(lb), options);

save('results_16_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_median_sigma_min')

%% Plot pareto front
load('results_16_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere.mat')
figure
scatter(-fval(:,1),-fval(:,2));
hold on
% scatter(0.012237, 1.62334e-7, '*')
grid on
xlabel('reciprocal condition number')
ylabel('min sigma_m')
title('Pareto front of the objectives')
legend('Pareto front', 'Grid result')

%%
[skew([1;0;0])*[0;0;1], skew([0;1;0])*[0;0;1], skew([0;0;1])*[0;0;1]]
skew([0;0;1])
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
    average_svd_set = [];

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

        % Exclude the 0 sigma
        sigma(end) = [];
        average_svd_set = [average_svd_set; mean(sigma)];
    end
    
    % Minimize the negative of the min -> maximize the min
    obj = [-min(reciprocal_condition_number_set), -median(min_svd_set)];
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


