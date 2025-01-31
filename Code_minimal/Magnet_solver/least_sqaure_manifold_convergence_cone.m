close all
clear 
clc

%% Constant
scale = 0.1;
sph_dia = 3.175e-3*[0; 0; 1]/scale;
Volumn = 24*(4/3)*pi*(norm(sph_dia)/2)^3;
B_r = 1.32;

type = "3D";

plot_progress = false;
iteration = 1000;
%% Magnet config in spherical coordinate
phi_sphere = linspace(0, 2*pi, 9);
phi_sphere(end) = [];
theta_sphere = linspace(0, pi, 7);
theta_sphere(1) = [];
theta_sphere(end) = [];
theta = linspace(-pi/2+deg2rad(15),pi/2-deg2rad(15),5);
phi = linspace(-pi/2+deg2rad(15), pi/2-deg2rad(15),5);

% 0+deg2rad(30),pi-deg2rad(30),5
% pi/2+deg2rad(30),3*pi/2-deg2rad(30),5
% pi+deg2rad(30),4*pi/2-deg2rad(30),5
% -pi/4+deg2rad(30),3*pi/4-deg2rad(30),5
% pi/4+deg2rad(30),5*pi/4-deg2rad(30),5
% 3*pi/4+deg2rad(30),7*pi/4-deg2rad(30),5
% 5*pi/4+deg2rad(30),9*pi/4-deg2rad(30),5

% theta = linspace(-pi,pi,9);
% theta(end) = [];
% phi = linspace(-pi/2,pi/2,9);
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

r_array = linspace(0,5,6)*0.01; % 5*sqrt(3) cm
r_array(1) = [];

for i = 1:length(r_array)
    for j = 1:size(sphere, 2)
        phi = sphere(1, j);
        theta = sphere(2, j);
        r = r_array(i);
        x = r*sin(theta)*cos(phi)/scale;
        y = r*sin(theta)*sin(phi)/scale;
        z = (r*cos(theta)+0.15)/scale;
    
        points = [points, [x;y;z;sphere(3,j);sphere(4,j);sphere(5,j)]];
    end
end

magnet_conf = points;

% pattern = [0;0;1.8];
% tolerance = 1e-2;
% matchingColumns = all(abs(magnet_conf(1:3, :) - pattern) < tolerance, 1);
% magnet_conf = magnet_conf(:,matchingColumns);

 %% Construct the ground truth magnetic field
% Sensor positions
% sens_pos_collection = [0.05, 0,0,0,0.05,0,-0.05,0,0,0,-0.05,0]/scale;
% sens_pos_collection = [0.10, 0,0,0,0.10,0,-0.10,0,0,0,-0.10,0]/scale;
% sens_pos_collection = [0.20, 0,0,0,0.20,0,-0.20,0,0,0,-0.20,0]/scale;
sens_pos_collection = [0.30, 0,0,0,0.30,0,-0.30,0,0,0,-0.30,0]/scale;
% sens_pos_collection = [1, 0,0,-1,0,0,0,1,0,0,-1,0]/scale;

sens_pos_collection = [sens_pos_collection, 4];


p_star_collection = {};
R_star_collection = {};
B_star_collection = {};


for magnet_conf_index = 1:size(magnet_conf,2)
    one_magnet_conf = magnet_conf(:, magnet_conf_index);
    p_star = one_magnet_conf(1:3);
    
    theta_star = one_magnet_conf(4);
    phi_star = one_magnet_conf(5);
    psi_star = one_magnet_conf(6);

    Rx_1 = [1 0 0;                  % rotate about x-axis
    0 cos(theta_star) -sin(theta_star);
    0 sin(theta_star) cos(theta_star)];
    
    Ry = [cos(phi_star) 0 sin(phi_star);    % rotate about y-axis
    0 1 0;
    -sin(phi_star) 0 cos(phi_star)];

    Rx_2 = [1 0 0;                  % rotate about x-axis
    0 cos(psi_star) -sin(psi_star);
    0 sin(psi_star) cos(psi_star)];

    % R_star = Rx_1*Ry*Rx_2;
    R_star = Ry*Rx_1*Rx_2;

    B_star = mag_field_vector(sens_pos_collection, p_star, B_r, Volumn, R_star);

    p_star_collection{end+1} = p_star;
    R_star_collection{end+1} = R_star;
    B_star_collection{end+1} = B_star;
end

B_star_corrupted_collection = {};
% Add sensor noise
for i = 1:length(B_star_collection)
    B_star_corrupted_collection{end+1} = B_star_collection{i} + 1e-7*randn(size(B_star_collection{i}));
end

% Run the Gauss-Newton Algorithm
stepsize_p = 0.08:0.005:0.08;
stepsize_R = pi/2:pi/10:pi/2;

result = [];

for i = 1:length(stepsize_p)
    for j = 1:length(stepsize_R)
        [p_k_collection, R_k_collection, B_k_collection, ...
         residual_B_collection, residual_p_collection, residual_R_collection, residual_R_angle_collection, ...
         num_not_correct_p, num_not_correct_R, mean_residual_p, mean_residual_R, elasped_time] = ...
        Gauss_Newton(stepsize_p(i), stepsize_R(j), ...
                     sens_pos_collection, magnet_conf, ...
                     p_star_collection, R_star_collection, B_star_collection, ...
                     100, plot_progress, scale, B_r, Volumn, type);

        % Count those are NaNs
        nan_columns = all(isnan([p_k_collection{:}]), 1);
        sum_NaN = sum(nan_columns);

        result = [result, [stepsize_p(i);stepsize_R(j);num_not_correct_p+sum_NaN;num_not_correct_R+sum_NaN; mean_residual_p; mean_residual_R; elasped_time]];
        fprintf('"stepsize_p": %.2f\n', stepsize_p(i));
        fprintf('"stepsize_R": %.2f\n', stepsize_R(j));
        fprintf('"Wrong Position": %d\n', num_not_correct_p + sum_NaN);
        fprintf('"Wrong Orientation": %d\n', num_not_correct_R + sum_NaN);
        fprintf('"Mean error p": %.2f\n', mean_residual_p);
        fprintf('"Mean error R": %.2f\n', mean_residual_R);
        fprintf('"Time used": %.2f seconds\n', elasped_time);

    end
end



% Gauss Newton Algorithm
function [p_k_collection, R_k_collection, B_k_collection, ...
          residual_B_collection, residual_p_collection, residual_R_collection, residual_R_angle_collection, ...
          num_not_correct_p, num_not_correct_R, mean_residual_p, mean_residual_R, elasped_time] = ...
         Gauss_Newton(stepsize_p, stepsize_R, ...
         sens_pos_collection, magnet_conf, ...
         p_star_collection, R_star_collection, B_star_collection, ...
         iteration, plot_progress, scale, B_r, Volumn, type)

    if plot_progress
        % Define the vertices of the cube
        vertices = [
            -0.05, -0.05, 0.10;  % Vertex 1
            -0.05,  0.05, 0.10;  % Vertex 2
             0.05, -0.05, 0.10;  % Vertex 3
             0.05,  0.05, 0.10;  % Vertex 4
            -0.05, -0.05, 0.20;  % Vertex 5
            -0.05,  0.05, 0.20;  % Vertex 6
             0.05, -0.05, 0.20;  % Vertex 7
             0.05,  0.05, 0.20;  % Vertex 8
        ]/scale;
    
        % Define the edges of the cube by connecting the vertices
        edges = [
            1, 2;
            1, 3;
            2, 4;
            3, 4;
            5, 6;
            5, 7;
            6, 8;
            7, 8;
            1, 5;
            2, 6;
            3, 7;
            4, 8;
        ];
    
        figure
        hold on
    
        % Plot the edges of the cube
        for i = 1:size(edges, 1)
            plot3([vertices(edges(i, 1), 1), vertices(edges(i, 2), 1)], ...
                  [vertices(edges(i, 1), 2), vertices(edges(i, 2), 2)], ...
                  [vertices(edges(i, 1), 3), vertices(edges(i, 2), 3)], 'k', 'LineWidth', 2);
        end
    
    
        axis equal
    end
    
    p_k_collection = {};
    R_k_collection = {};
    B_k_collection = {};
    residual_B_collection = {};
    
    % Initial guess
    multi_start_p = {[-0.0;0.0;0.15]/scale};
    multi_start_R = {rotx(0)};

    % multi_start_p = {[-0.63629;-1.03374;0.22846]/scale};
    % multi_start_R = {[0.9952 -0.0342 -0.0922; -0.0581 0.5507 -0.8326; 0.0793 0.8339 0.5461]};

    lb_p = [-0.05, -0.05, 0.10]/scale;
    ub_p = [0.05, 0.05, 0.20]/scale;
    
    tic
    % Solve for each magnet configuration
    for magnet_conf_index = 1:size(magnet_conf,2)
        
        % Plot the ground truth config
        if plot_progress
            one_magnet_conf = magnet_conf(:, magnet_conf_index);
            p_star = one_magnet_conf(1:3);
            
            theta_star = one_magnet_conf(4);
            phi_star = one_magnet_conf(5);
            psi_star = one_magnet_conf(6);
        
            Rx_1 = [1 0 0;                  % rotate about x-axis
            0 cos(theta_star) -sin(theta_star);
            0 sin(theta_star) cos(theta_star)];
            
            Ry = [cos(phi_star) 0 sin(phi_star);    % rotate about y-axis
            0 1 0;
            -sin(phi_star) 0 cos(phi_star)];
        
            Rx_2 = [1 0 0;                  % rotate about x-axis
            0 cos(psi_star) -sin(psi_star);
            0 sin(psi_star) cos(psi_star)];
        
            R_star = Rx_1*Ry*Rx_2;
    
            u = [0;0;1];
            direction = R_star*u;
        
            arrow_length = 0.5;
            direction = direction*arrow_length;
            
            quiver3(p_star(1), p_star(2), p_star(3), direction(1), direction(2), direction(3), ...
                'LineWidth', 2, 'MaxHeadSize', 100, 'Color', 'r');   
    
        end
    
        residual_multistart_collection = {};
        p_k_multistart_collection = {};
        R_k_multistart_collection = {};
        B_k_multistart_collection = {};
        
        % Extract the ground truth magnetic field
        B_star = B_star_collection{magnet_conf_index};
    
        % Constraint the step size such that the local linearization is valid
        max_change_norm_p = stepsize_p;
        max_change_norm_R = stepsize_R;
        
        % Start from multiple initial guess
        for initial_p_index = 1:length(multi_start_p)
            for initial_R_index = 1:length(multi_start_R)
                initial_guess_p = multi_start_p{initial_p_index};
                initial_guess_R = multi_start_R{initial_R_index};
        
                for iter = 1:iteration  
                    % iter
                    % Exit condition 1: estimated field is close enough to the
                    % ground truth
                    if iter == 1
                        B_estimation = mag_field_vector(sens_pos_collection, initial_guess_p, B_r, Volumn, initial_guess_R);
                        residual_initial = norm(B_estimation-B_star);
                        if norm(B_estimation-B_star) <= 1e-20
                            break
                        end
                    else
                        B_estimation = mag_field_vector(sens_pos_collection, p_k, B_r, Volumn, R_k);
                        if norm(B_estimation-B_star) <= 1e-20
                            % i
                            break
                        end
                    end
            
            
                    if iter == 1
                        if plot_progress
                            position = initial_guess_p;
                            rotation = initial_guess_R;
                            
                            u = [0;0;1];
                            direction = rotation*u;
                        
                            arrow_length = 0.5;
                            direction = direction*arrow_length;
                            
                            quiver3(position(1), position(2), position(3), direction(1), direction(2), direction(3), ...
                                'LineWidth', 2, 'MaxHeadSize', 100, 'Color', 'b');         
                        end
    
                        update = step(B_star, sens_pos_collection, initial_guess_p, B_r, Volumn, initial_guess_R, type);
                        delta_p = update(1:3);
                        delta_R = update(4:6);
    
                                            
                        % Exit condition 2: step size too small
                        if norm(delta_p) <= 1e-12 && norm(delta_R) <= 1e-12
                            % i
                            break
                        end
    
                        % Prevent numerical issue
                        if norm(delta_p) <= 1e-12
                            delta_p = [0;0;0];
                        end
                        if norm(delta_R) <= 1e-12
                            delta_R = [0;0;0];
                        end
    
    
                        % Constraint the step size such that the norms of delta_p_scaled and
                        % delta_R_scaled are <= max_change_norm
                        delta_p_norm = norm(delta_p);
                        delta_R_norm = norm(delta_R);

                        % if delta_p_norm >= max_change_norm_p || delta_R_norm >= max_change_norm_R
                        %     ratio_p = delta_p_norm/max_change_norm_p;
                        %     ratio_R = delta_R_norm/max_change_norm_R;
                        % 
                        %     ratio = max([ratio_p, ratio_R]);
                        %     c_p = 1/ratio;
                        %     c_R = 1/ratio;
                        % else
                        %    c_p = 1;
                        %    c_R = 1;
                        % end

                        % if delta_p_norm >= max_change_norm_p
                        %    c_p = max_change_norm_p/delta_p_norm; 
                        % else
                        %    c_p = 1;
                        % end
                        % 
                        % if delta_R_norm >= max_change_norm_R
                        %     c_R = max_change_norm_R/delta_R_norm;
                        % else
                        %     c_R = 1;
                        % end

                        c_p = 1;
                        c_R = 1;

                        
                        p_k = initial_guess_p + c_p *delta_p;
                        R_k = initial_guess_R*expm(skew(c_R*delta_R));
    
                        % % Project the position back to the limits
                        % for j = 1:length(p_k)
                        %     p_k(j) = max(lb_p(j), min(ub_p(j), p_k(j)));
                        % end
                    else
                        if plot_progress
                            position = p_k;
                            rotation = R_k;
                            
                            u = [0;0;1];
                            direction = rotation*u;
                        
                            arrow_length = 0.5;
                            direction = direction*arrow_length;
                            
                            quiver3(position(1), position(2), position(3), direction(1), direction(2), direction(3), ...
                                'LineWidth', 2, 'MaxHeadSize', 100, 'Color', 'b');  
                        end
    
                        update = step(B_star, sens_pos_collection, p_k, B_r, Volumn, R_k, type);
                        delta_p = update(1:3);
                        delta_R = update(4:6);
                        
        
                        % Exit condition 2: step size too small
                        if norm(delta_p) <= 1e-12 && norm(delta_R) <= 1e-12
                            break
                        end
    
                        % Prevent numerical issue
                        if norm(delta_p) <= 1e-12
                            delta_p = [0;0;0];
                        end
                        if norm(delta_R) <= 1e-12
                            delta_R = [0;0;0];
                        end
    
    
                        % Constraint the step size such that the norms of delta_p_scaled and
                        % delta_R_scaled are <= max_change_norm
                        delta_p_norm = norm(delta_p);
                        delta_R_norm = norm(delta_R);

                        % if delta_p_norm >= max_change_norm_p || delta_R_norm >= max_change_norm_R
                        %     ratio_p = delta_p_norm/max_change_norm_p;
                        %     ratio_R = delta_R_norm/max_change_norm_R;
                        % 
                        %     ratio = max([ratio_p, ratio_R]);
                        %     c_p = 1/ratio;
                        %     c_R = 1/ratio;
                        % else
                        %    c_p = 1;
                        %    c_R = 1;
                        % end

                        % if delta_p_norm >= max_change_norm_p
                        %    c_p = max_change_norm_p/delta_p_norm; 
                        % else
                        %    c_p = 1;
                        % end
                        % 
                        % if delta_R_norm >= max_change_norm_R
                        %     c_R = max_change_norm_R/delta_R_norm;
                        % else
                        %     c_R = 1;
                        % end

                        c_p = 1;
                        c_R = 1;

                        p_k = p_k + c_p *delta_p;
                        R_k = R_k*expm(skew(c_R*delta_R));
    
                        % % Project the position back to the limits
                        % for j = 1:length(p_k)
                        %     p_k(j) = max(lb_p(j), min(ub_p(j), p_k(j)));
                        % end
                        
                    end
                end
                
                % Calculate the magnetic field corresponding to the estimated
                % magnet config
                B_k = mag_field_vector(sens_pos_collection, p_k, B_r, Volumn, R_k);
    
                % Put the results into the collection
                residual_multistart_collection{end+1} = norm(B_k-B_star);
                p_k_multistart_collection{end+1} = p_k;
                R_k_multistart_collection{end+1} = R_k;
                B_k_multistart_collection{end+1} = B_k;
        
            end
        end
        
        % Find the one that possibly converges
        residual_multistart_collection_array = cell2mat(residual_multistart_collection);
        [smallest_residual, index] = min(residual_multistart_collection_array);
    
        % Put the possibly converged result into the collection
        p_k_collection{end+1} = p_k_multistart_collection{index};
        R_k_collection{end+1} = R_k_multistart_collection{index};
        B_k_collection{end+1} = B_k_multistart_collection{index};
        residual_B_collection{end+1} = residual_multistart_collection{index};
    end
    elasped_time = toc;
    
    % Calculate the errors
    residual_p_collection = [];
    residual_R_collection = [];
    residual_R_angle_collection = [];
    for config_num = 1:size(magnet_conf,2)
        residual_p_collection = [residual_p_collection, norm(p_k_collection{config_num}-p_star_collection{config_num})];
    
        R_k = R_k_collection{config_num};
        R_star = R_star_collection{config_num};
        residual_R_collection = [residual_R_collection, norm(R_k(:,3)-R_star(:,3))];
        cos_theta = dot(R_k(:,3),R_star(:,3));
        theta = acos(cos_theta);
        residual_R_angle_collection = [residual_R_angle_collection, rad2deg(theta)];
    end
    
    num_not_correct_p = sum(residual_p_collection>0.1);
    num_not_correct_R = sum(residual_R_angle_collection>5);

    mean_residual_p = mean(residual_p_collection);
    mean_residual_R = mean(residual_R_collection);

end

%% Incremental step in least square algorithm
function update = step(measurement, sens_pos_collection, p_k, B_r, Volumn, R_k, type)

    % Estimated field for current magnet's pose's estimation
    estimation = mag_field_vector(sens_pos_collection, p_k, B_r, Volumn, R_k);
    
    % Error vector
    error = (measurement - estimation);
    
    % Construct sensor config, position
    sens_num = sens_pos_collection(end);
    sens_pos_collection(end) = [];
    sens_pos_collection = reshape(sens_pos_collection, 3, []);
    
    % Construct sensor config, orientation
    % Default orientation
    default_or = [1;0;0;0];

    % Orientation for all sensors
    sens_or = repmat(default_or, 1, sens_num);
    
    % Construct the Jacobian form each sensor
    J = [];
    for i = 1:sens_num
        J = [J; J_analytical_sensor_B_to_world_pR(sens_pos_collection(:,i), sens_or(:,i), ...
            p_k, B_r, Volumn, R_k, type)];
    end
    
    % Incremental step is [delta_p; delta_R] = pinv(J)*(-error)
    update = pinv(J)*(error);

    % if norm(estimation+measurement)<=1e-15
    %     delta_R = [pi;0;0];
    %     update(4:6) = delta_R;
    % end
end

function residual = residual_func(measurement, sens_pos_collection, p_k, B_r, Volumn, R_k)
    % Estimated field for current magnet's pose's estimation
    estimation = mag_field_vector(sens_pos_collection, p_k, B_r, Volumn, R_k);
    
    % Error vector
    residual = (measurement - estimation);
end

%% Quaternion to rotation matrix
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

%% Skew symmetric matrix
% Convert a vector to a skew symmetric matrix
function v_hat = skew(v)
    v_hat = [0 -v(3) v(2); 
             v(3) 0 -v(1);
             -v(2) v(1) 0];
end

%% Magnetic field vector for all sensors
function B_vector = mag_field_vector(sens_pos_collection, magnet_pos, B_r, Volumn, R_star)
    B_vector = [];
    sens_num = sens_pos_collection(end);
    sens_pos_collection(end) = [];
    sens_pos_collection  = reshape(sens_pos_collection, 3, []);
    for sens_index = 1:sens_num
        sens_pos = sens_pos_collection(:,sens_index);
        B_vector = [B_vector; mag_field(sens_pos, magnet_pos, B_r, Volumn, R_star)];
    end
end

%% Magnetic field
function B = mag_field(sens_pos, magnet_pos, B_r, Volumn, R_star)
    r = sens_pos-magnet_pos;
    r_hat = r/norm(r);
    B = (B_r*Volumn)/(4*pi*(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*R_star*[0;0;1];
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