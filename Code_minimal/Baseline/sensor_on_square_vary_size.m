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

%% 2 by 2 grid
side_length = 0.02:0.02:0.5;

rcond_all_sens_conf = [];
min_svd_all_sens_conf = [];
median_min_svd_all_sens_conf = [];
mean_min_svd_all_sens_conf = [];
max_svd_all_sens_conf = [];
mean_svd_all_sens_conf = [];
for i = 1:length(side_length)
    sens_conf = [side_length(i)/scale side_length(i)/scale 0 1 0 0 0 ...
                 side_length(i)/scale -side_length(i)/scale 0 1 0 0 0 ...
                 -side_length(i)/scale side_length(i)/scale 0 1 0 0 0 ...
                 -side_length(i)/scale -side_length(i)/scale 0 1 0 0 0 ...
                 4];

    sens_num = sens_conf(end);
    sens_conf(end) = [];
    sens_conf = reshape(sens_conf, 7, []);
    
    % Extract position and orientation
    sens_pos = sens_conf(1:3,:);
    
    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;

    rcond_one_sens_conf = [];
    min_svd_one_sens_conf = [];
    max_svd_one_sens_conf = [];
    mean_svd_one_sens_conf = [];

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

        % Exclude the 0 sigma
        sigma_scaled(end) = [];
        mean_svd_one_sens_conf = [mean_svd_one_sens_conf; mean(sigma_scaled)];
    end
    
    % Put the min value in the workspace into the list
    rcond_all_sens_conf = [rcond_all_sens_conf; min(rcond_one_sens_conf)];
    min_svd_all_sens_conf = [min_svd_all_sens_conf; min(min_svd_one_sens_conf)];
    median_min_svd_all_sens_conf = [median_min_svd_all_sens_conf; median(min_svd_one_sens_conf)];
    mean_min_svd_all_sens_conf = [mean_min_svd_all_sens_conf; mean(min_svd_one_sens_conf)];
    max_svd_all_sens_conf = [max_svd_all_sens_conf; min(max_svd_one_sens_conf)];
    mean_svd_all_sens_conf = [mean_svd_all_sens_conf; mean(mean_svd_one_sens_conf)];
end

%% Plots
figure(1)
plot(2*side_length, rcond_all_sens_conf, "linewidth", 2);
% ylim([0 1])
xlabel("Side length of the square (m)")
ylabel("$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$", "Interpreter", "latex")
grid on

% figure(2)
% plot(2*side_length, min_svd_all_sens_conf, "linewidth", 2);
% % ylim([0 2.5e-7])
% xlabel("Side length of the square (m)")
% ylabel("Minimum of the \sigma_{min} in the workspace")
% grid on

% figure(3)
% plot(2*side_length, mean_svd_all_sens_conf, "linewidth", 2);
% % ylim([0 2.5e-7])
% xlabel("Side length of the square (m)")
% ylabel("Mean sigma in the workspace")
% grid on

figure(4)
plot(2*side_length, median_min_svd_all_sens_conf, "linewidth", 2);
% ylim([0 2.5e-7])
xlabel("Side length of the square (m)", "Interpreter", "latex")
ylabel("$median_{x \in \mathcal{X}} \sigma_{min}(x)$", "Interpreter", "latex")
grid on

figure(5)
plot(2*side_length, mean_min_svd_all_sens_conf, "linewidth", 2);
% ylim([0 2.5e-7])
xlabel("Side length of the square (m)", "Interpreter", "latex")
ylabel("$mean_{x \in \mathcal{X}} \sigma_{min}(x)$", "Interpreter", "latex")
grid on

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
