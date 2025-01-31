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

%% Read the experiment data
% Folder path
folder = fullfile(pwd, "data_10cm");

% .mat files
matFiles = dir(fullfile(folder, '*.mat'));

B_star = [];
magnet_conf = [];
B_star_multiple_time_stamps = {};

% Go through the data
for i = 1:length(matFiles)
    % Get the full path of the .mat file
    matFilePath = fullfile(folder, matFiles(i).name);

    % Load the .mat file
    data = load(matFilePath);
    
    % Process the loaded data (replace this with your actual processing code)
    disp(['Processing file: ', matFiles(i).name]);

    % Extract the filename without the extension, for the magnet
    % configuration
    [~, fileName, ~] = fileparts(matFiles(i).name);
    
    % Split the filename by the underscore "_"
    parts = strsplit(fileName, '_');
    
    % Extract the first part and the second part (if they exist)
    if length(parts) >= 2
        firstPart = parts{1};  % Part before the first underscore
        secondPart = parts{2}; % Part between the first and second underscores
    elseif length(parts) == 1
        firstPart = parts{1};  % If no underscore, entire name is the first part
        secondPart = '';       % No second part available
    else
        firstPart = '';
        secondPart = '';
    end
    
    % Display the extracted parts
    disp(['Position: ', firstPart]);
    disp(['Y orientation: ', secondPart]);

    if firstPart == "0010"
        p_m = [0;0;0.1];
    elseif  firstPart == "0015"
        p_m = [0;0;0.15];
    elseif firstPart == "0020"
        p_m = [0;0;0.2];
    elseif firstPart == "0n515"
        p_m = [0;-0.05;0.15];
    elseif firstPart == "n5015"
        p_m = [-0.05;0;0.15];
    elseif firstPart == "5015"
        p_m = [0.05;0;0.15];
    elseif firstPart == "0515"
        p_m = [0;0.05;0.15];
    end

    if secondPart == "y0"
        R_y = deg2rad(0);
    elseif secondPart == "y45"
        R_y = deg2rad(45);
    elseif secondPart == "y90"
        R_y = deg2rad(90);
    elseif secondPart == "y135"
        R_y = deg2rad(135);
    end

    
    % Go through each saved field
    for j = 1:length(data.save_element)
        B_star = [B_star, mean(data.user_data.Field_filt(:, data.save_element(j)-11: data.save_element(j)), 2)];
        B_star_multiple_time_stamps{end+1} = data.user_data.Field_filt(:, data.save_element(j)-11: data.save_element(j));

        if j == 1
            R_x = deg2rad(0);
        elseif j == 2
            R_x = deg2rad(45);
        elseif j == 3
            R_x = deg2rad(90);
        elseif j == 4
            R_x = deg2rad(135);
        elseif j == 5
            R_x = deg2rad(180);
        elseif j == 6
            R_x = deg2rad(225);
        elseif j == 7
            R_x = deg2rad(270);
        elseif j == 8
            R_x = deg2rad(315);
        end

        R = [R_y; R_x; 0];
        
        magnet_conf = [magnet_conf, [p_m;R]];

    end
end
%% Construct the ground truth magnetic field
% Sensor positions
% sens_pos_collection = [0.05, 0,0,0,0.05,0,-0.05,0,0,0,-0.05,0]/scale;
% sens_pos_collection = [0.0459,-0.0016,0.0017, ...
%                        -0.0015,0.0503,0.0028, ...
%                        -0.0518,-0.0009,0.0013, ...
%                        -0.0004,-0.0500,0.0007]/scale; % fminimax
% sens_pos_collection = [0.0463,-0.0012,0.0015, ...
%                        -0.0039,0.0507,0.0037, ...
%                        -0.0537,-0.0015,0.0025, ...
%                        -0.0040,-0.0515,0.0011]/scale; % lsqnonlin

% sens_pos_collection = [0.10, 0,0,0,0.10,0,-0.10,0,0,0,-0.10,0]/scale;
sens_pos_collection = [0.0972,-0.0013,0.0002, ...
                       -0.0039,0.1015,0.0030, ...
                       -0.1054,-0.0017,0.0029, ...
                       -0.0038,-0.1017,-0.0003]/scale;


% sens_pos_collection = [0.1120, 0,0,-0.1120,0,0,0,0.1120,0,0,-0.1120,0]/scale;
% sens_pos_collection = [0.2150, 0,0,-0.2150,0,0,0,0.2150,0,0,-0.2150,0]/scale;
% sens_pos_collection = [0.3550, 0,0,-0.3550,0,0,0,0.3550,0,0,-0.3550,0]/scale;
% sens_pos_collection = [1, 0,0,-1,0,0,0,1,0,0,-1,0]/scale;



sens_pos_collection = [sens_pos_collection, 4];


p_star_collection = {};
R_star_collection = {};
B_star_collection = {};
B_star_multiple_time_stamps_collection = {};

for magnet_conf_index = 1:size(magnet_conf,2)
    one_magnet_conf = magnet_conf(:, magnet_conf_index);
    p_star = one_magnet_conf(1:3);
   
    phi_star = one_magnet_conf(4);
    theta_star = one_magnet_conf(5);

    Ry = [cos(phi_star) 0 sin(phi_star);    % rotate about y-axis
    0 1 0;
    -sin(phi_star) 0 cos(phi_star)];

    Rx = [1 0 0;                  % rotate about x-axis
    0 cos(theta_star) -sin(theta_star);
    0 sin(theta_star) cos(theta_star)];

    R_star = Ry*Rx;

    p_star_collection{end+1} = p_star/scale;
    R_star_collection{end+1} = R_star;
    B_star_collection{end+1} = B_star(:, magnet_conf_index)*1e-6;
    B_star_multiple_time_stamps_collection{end+1} = B_star_multiple_time_stamps{magnet_conf_index}*1e-6;
end

% p_star_collection = {p_star_collection{2}};
% R_star_collection = {R_star_collection{2}};
% B_star_collection = {B_star_collection{2}};
% magnet_conf = magnet_conf(:,2);

%% Multiple measurement for one configuration
mean_residual_p = [];
mean_residual_R = [];
for i = 1:length(B_star_multiple_time_stamps_collection)
    B_collection_one_config = cell(1, size(B_star_multiple_time_stamps_collection{i}, 2));
    % Loop through each column and store it in the cell array
    for j = 1:size(B_star_multiple_time_stamps_collection{i}, 2)
        B_collection_one_config{j} = B_star_multiple_time_stamps_collection{i}(:, j);  % Store each column as a separate element in the cell array
    end

    p_star_collection_one_config = cell(1, size(B_star_multiple_time_stamps_collection{i}, 2));
    for j = 1:size(B_star_multiple_time_stamps_collection{i}, 2)
        p_star_collection_one_config{j} = p_star_collection{i};
    end

    R_star_collection_one_config = cell(1, size(B_star_multiple_time_stamps_collection{i}, 2));
    for j = 1:size(B_star_multiple_time_stamps_collection{i}, 2)
        R_star_collection_one_config{j} = R_star_collection{i};
    end


    magnet_conf_one_config = repmat(magnet_conf(:,i), 1, size(B_collection_one_config,2));
    
    stepsize_p = 0;
    stepsize_R = 0;
    [p_k_collection, R_k_collection, B_k_collection, ...
         residual_B_collection, residual_p_collection, residual_R_collection, residual_R_angle_collection, ...
         num_not_correct_p, num_not_correct_R, elasped_time] = ...
        Gauss_Newton(stepsize_p, stepsize_R, ...
                     sens_pos_collection, magnet_conf_one_config, ...
                     p_star_collection_one_config, R_star_collection_one_config, B_collection_one_config, ...
                     100, plot_progress, scale, B_r, Volumn, type);

    mean_residual_p = [mean_residual_p, mean(residual_p_collection)];
    mean_residual_R = [mean_residual_R, mean(residual_R_angle_collection)];
end
mean(mean_residual_p)
mean(mean_residual_R)

%% Run the Gauss-Newton Algorithm
stepsize_p = 0.08:0.005:0.08;
stepsize_R = pi/2:pi/10:pi/2;

result = [];

for i = 1:length(stepsize_p)
    for j = 1:length(stepsize_R)
        [p_k_collection, R_k_collection, B_k_collection, ...
         residual_B_collection, residual_p_collection, residual_R_collection, residual_R_angle_collection, ...
         num_not_correct_p, num_not_correct_R, elasped_time] = ...
        Gauss_Newton(stepsize_p(i), stepsize_R(j), ...
                     sens_pos_collection, magnet_conf, ...
                     p_star_collection, R_star_collection, B_star_collection, ...
                     100, plot_progress, scale, B_r, Volumn, type);

        % Count those are NaNs
        nan_columns = all(isnan([p_k_collection{:}]), 1);
        sum_NaN = sum(nan_columns);

        result = [result, [stepsize_p(i);stepsize_R(j);num_not_correct_p+sum_NaN;num_not_correct_R+sum_NaN;elasped_time]]
    end
end



% Gauss Newton Algorithm
function [p_k_collection, R_k_collection, B_k_collection, ...
          residual_B_collection, residual_p_collection, residual_R_collection, residual_R_angle_collection, ...
          num_not_correct_p, num_not_correct_R, elasped_time] = ...
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
    
    lb_p = [-0.05, -0.05, 0.10]/scale;
    ub_p = [0.05, 0.05, 0.20]/scale;
    
    tic
    % Solve for each magnet configuration
    for magnet_conf_index = 1:size(magnet_conf,2)

        % Initial guess
        multi_start_p = {p_star_collection{magnet_conf_index}};
        multi_start_R = {R_star_collection{magnet_conf_index}};
        
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

end

%% Projected error
sens_num = sens_pos_collection(end);
sens_pos_collection(end) = [];
sens_pos_collection = reshape(sens_pos_collection, 3, sens_num);

% Default orientation
default_or = [1;0;0;0];

% Orientation for all sensors
sens_or_unitary = repmat(default_or, 1, sens_num);

singular_value_with_projected_error = {};
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
        J_scaled = [J_scaled;J_analytical_sensor_B_to_world_pR(sens_pos_collection(:,j), sens_or_unitary(:,j), ...
            magnet_pos, B_r, Volumn, R_star, type)];
    end

    [U,S,V] = svd(J_scaled);

    num_dof = 5;
    
    v1 = V(:,1);
    v2 = V(:,2);
    v3 = V(:,3);
    v4 = V(:,4);
    v5 = V(:,5);

    p_k = p_k_collection{magnet_num};
    R_k = R_k_collection{magnet_num};
    R_0 = R_star_collection{magnet_num};
    p_0 = p_star_collection{magnet_num};
    R_rel = R_0.'*R_k;
    theta = acos((trace(R_rel)-1)/2);
    w_hat =(1/sin(theta))*[R_rel(3,2)-R_rel(2,3);R_rel(1,3)-R_rel(3,1);R_rel(2,1)-R_rel(1,2)];
    w=theta*w_hat;
    delta_p = p_k - p_0;
    error = [delta_p; w];
    norm(error);
    
    error_project_1 = dot(error,v1);
    error_project_2 = dot(error,v2);
    error_project_3 = dot(error,v3);
    error_project_4 = dot(error,v4);
    error_project_5 = dot(error,v5);

    singular_value_with_projected_error{end+1} = abs([S(1,1), S(2,2), S(3,3), S(4,4), S(5,5);
                                                  S(1,1)/S(1,1), S(2,2)/S(1,1), S(3,3)/S(1,1), S(4,4)/S(1,1), S(5,5)/S(1,1);
                                                  error_project_1, error_project_2, error_project_3, error_project_4, error_project_5;
                                                  error_project_1/norm(error), error_project_2/norm(error), error_project_3/norm(error), error_project_4/norm(error), error_project_5/norm(error)]);

end

% Initialize a matrix to store the sum of all 4x5 arrays
sumArray = zeros(4, 5);

% Loop over each cell and accumulate the sum of the 4x5 matrices
for i = 1:size(singular_value_with_projected_error, 2)
    sumArray = sumArray + singular_value_with_projected_error{i};
end

% Calculate the mean by dividing by the number of arrays (224)
meanArray = sumArray / size(singular_value_with_projected_error, 2);
%% Process the results
position_0010 = [0; 0; 0.1];
position_0015 = [0; 0; 0.15];
position_0020 = [0; 0; 0.20];
position_5015 = [0.05; 0; 0.15];
position_0515 = [0; 0.05; 0.15];
position_n5015 = [-0.05; 0; 0.15];
position_0n515 = [0; -0.05; 0.15];

columns_0010 = all(magnet_conf(1:3, :) == position_0010, 1);
columns_0015 = all(magnet_conf(1:3, :) == position_0015, 1);
columns_0020 = all(magnet_conf(1:3, :) == position_0020, 1);
columns_5015 = all(magnet_conf(1:3, :) == position_5015, 1);
columns_0515 = all(magnet_conf(1:3, :) == position_0515, 1);
columns_n5015 = all(magnet_conf(1:3, :) == position_n5015, 1);
columns_0n515 = all(magnet_conf(1:3, :) == position_0n515, 1);

residual_p_collection_0010 = residual_p_collection(:, columns_0010);
residual_R_collection_0010 = residual_R_collection(:, columns_0010);
residual_R_angle_collection_0010 = residual_R_angle_collection(:, columns_0010);


residual_p_collection_0015 = residual_p_collection(:, columns_0015);
residual_R_collection_0015 = residual_R_collection(:, columns_0015);
residual_R_angle_collection_0015 = residual_R_angle_collection(:, columns_0015);


residual_p_collection_0020 = residual_p_collection(:, columns_0020);
residual_R_collection_0020 = residual_R_collection(:, columns_0020);
residual_R_angle_collection_0020 = residual_R_angle_collection(:, columns_0020);


residual_p_collection_5015 = residual_p_collection(:, columns_5015);
residual_R_collection_5015 = residual_R_collection(:, columns_5015);
residual_R_angle_collection_5015 = residual_R_angle_collection(:, columns_5015);


residual_p_collection_0515 = residual_p_collection(:, columns_0515);
residual_R_collection_0515 = residual_R_collection(:, columns_0515);
residual_R_angle_collection_0515 = residual_R_angle_collection(:, columns_0515);


residual_p_collection_n5015 = residual_p_collection(:, columns_n5015);
residual_R_collection_n5015 = residual_R_collection(:, columns_n5015);
residual_R_angle_collection_n5015 = residual_R_angle_collection(:, columns_n5015);


residual_p_collection_0n515 = residual_p_collection(:, columns_0n515);
residual_R_collection_0n515 = residual_R_collection(:, columns_0n515);
residual_R_angle_collection_0n515 = residual_R_angle_collection(:, columns_0n515);






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