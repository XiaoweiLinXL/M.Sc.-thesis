%%
clear all
close all
clc

%% Plot the minimum condition number among the 5 optimized configurations Multiobj
max_min_rcond_set_multiobj = [];

for i = 2:10
    % Load the solution
    % Construct the file name
    fileName = sprintf('results_%d_3_axis_multiobj_2000gen_20000pop.mat', i);

    load(fileName);
    
    sens_conf = [sol];
    [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

    max_min_rcond = max(min_rcond);
    max_min_rcond_set_multiobj = [max_min_rcond_set_multiobj;max_min_rcond];

end

%% Plots Multiobj - objectives
sensor_number = [2;3;4;5;6;7;8;9;10];
plot(sensor_number*3, max_min_rcond_set_multiobj, '-x', 'LineWidth',2, 'MarkerSize',8)
ylim([0,0.1])
grid on

% Add text labels beside each data point
for i = 1:length(sensor_number)
    text(sensor_number(i)*3, max_min_rcond_set_multiobj(i), sprintf('%.4f', max_min_rcond_set_multiobj(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

hold on

%% Get the minimum rcond among the entire workspace for multiobj optimization
min_reciprocal_condition_number_set_multiobj = [];

for i = 2:10
    % Load the solution
    % Construct the file name
    fileName = sprintf('results_%d_3_axis_multiobj_2000gen_20000pop.mat', i);

    load(fileName);

    % Evaluate the sensor config
    sens_conf = [sol];
    [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

    obj_rcond

    % Get the sensor configuration that corresponds to the biggest minimum
    % condition number
    index = find(min_rcond == max(min_rcond));
    solution = sol(index,:);
    number_of_sensor = solution(end);
    solution(end) = [];
    sens_conf = reshape(solution, 7, []);
    
    % Extract position and orientation of the sensor configuration
    sens_pos = sens_conf(1:3,:);
    
    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Update the magnet configurations to a box and various orientations
    % Workspace as a box in m
    LB = [-0.05, -0.05, 0.2, deg2rad(-180), deg2rad(-90)];
    UB = [0.05, 0.05, 0.3, deg2rad(180), deg2rad(90)];
    
    % Number of samples within the box will be num_samp^3
    num_samp = 10;
    
    C = cell(1, 5);
    [C{:}] = ndgrid(linspace(0, 1, num_samp));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    magnet_conf = conf;

    % add a row of 0 for psi
    zero_row = zeros(1, size(magnet_conf,2));
    magnet_conf = [magnet_conf ; zero_row];
    
    % Get the rcond for all magnet configs in the workspace, get the
    % minimum of them
    reciprocal_condition_number_set = [];
    
    for magnet_num = 1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
        
        
        J = [];
        % Construct the Jacobian
        for sensor_num = 1:number_of_sensor
            J = [J; jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_X(sens_pos(:,sensor_num), sens_or_unitary(:,sensor_num), ...
                        magnet_pos, mu_norm, theta, phi, psi, type)];
        end
    
        sigma = svd(J);
        
        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);
    
        reciprocal_condition_number_set = [reciprocal_condition_number_set;reciprocal_condition_number];
    end
    
    min_reciprocal_condition_number_set_multiobj = [min_reciprocal_condition_number_set_multiobj; min(reciprocal_condition_number_set)];

end

%% Plots Multiobj - larger workspace
sensor_number = [2;3;4;5;6;7;8;9;10];
plot(sensor_number*3, min_reciprocal_condition_number_set_multiobj, '-bx', 'LineWidth',2, 'MarkerSize',8)
ylim([0,0.1])
grid on

% Add text labels beside each data point
for i = 1:length(sensor_number)
    text(sensor_number(i)*3, min_reciprocal_condition_number_set_multiobj(i), sprintf('%.4f', min_reciprocal_condition_number_set_multiobj(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

hold on


%% Plot the minimum condition number among the 5 optimized configurations Min rcond
max_min_rcond_set_minrcond = [];
for i = 2:9
    % Load the solution
    % Construct the file name
    fileName = sprintf('results_%d_3_axis_minrcond_2000gen_20000pop.mat', i);

    load(fileName);
    
    sens_conf = [sol];
    [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

    max_min_rcond = max(min_rcond);
    max_min_rcond_set_minrcond = [max_min_rcond_set_minrcond;max_min_rcond];

end

%% Plots min rcond - objectives
sensor_number = [2;3;4;5;6;7;8;9];
plot(sensor_number*3, max_min_rcond_set_minrcond, ':ro', 'LineWidth',2, 'MarkerSize',8)
ylim([0,0.1])
grid on


% Add text labels beside each data point
for i = 1:length(sensor_number)
    if sensor_number(i) == 7
        text(sensor_number(i)*3, max_min_rcond_set_minrcond(i)+0.005, sprintf('%.4f', max_min_rcond_set_minrcond(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    else
        text(sensor_number(i)*3, max_min_rcond_set_minrcond(i), sprintf('%.4f', max_min_rcond_set_minrcond(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end

hold on





%% Get the minimum rcond among the entire workspace for min rcond optimization
min_reciprocal_condition_number_set_minrcond = [];

for i = 2:9
    % Load the solution
    % Construct the file name
    fileName = sprintf('results_%d_3_axis_minrcond_2000gen_20000pop.mat', i);

    load(fileName);

    % Evaluate the sensor config
    sens_conf = [sol];
    [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

    % Get the sensor configuration that corresponds to the biggest minimum
    % condition number
    index = find(min_rcond == max(min_rcond));
    solution = sol(index,:);
    number_of_sensor = solution(end);
    solution(end) = [];
    sens_conf = reshape(solution, 7, []);
    
    % Extract position and orientation of the sensor configuration
    sens_pos = sens_conf(1:3,:);
    
    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Update the magnet configurations to a box and various orientations
    % Workspace as a box in m
    LB = [-0.05, -0.05, 0.2, deg2rad(-180), deg2rad(-90)];
    UB = [0.05, 0.05, 0.3, deg2rad(180), deg2rad(90)];
    
    % Number of samples within the box will be num_samp^3
    num_samp = 10;
    
    C = cell(1, 5);
    [C{:}] = ndgrid(linspace(0, 1, num_samp));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    magnet_conf = conf;

    % add a row of 0 for psi
    zero_row = zeros(1, size(magnet_conf,2));
    magnet_conf = [magnet_conf ; zero_row];
    
    % Get the rcond for all magnet configs in the workspace, get the
    % minimum of them
    reciprocal_condition_number_set = [];
    
    for magnet_num = 1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
        
        
        J = [];
        % Construct the Jacobian
        for sensor_num = 1:number_of_sensor
            J = [J; jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_X(sens_pos(:,sensor_num), sens_or_unitary(:,sensor_num), ...
                        magnet_pos, mu_norm, theta, phi, psi, type)];
        end
    
        sigma = svd(J);
        
        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);
    
        reciprocal_condition_number_set = [reciprocal_condition_number_set;reciprocal_condition_number];
    end
    
    min_reciprocal_condition_number_set_minrcond = [min_reciprocal_condition_number_set_minrcond; min(reciprocal_condition_number_set)];

end

%% Plots Min rcond - larger workspace
sensor_number = [2;3;4;5;6;7;8;9];
plot(sensor_number*3, min_reciprocal_condition_number_set_minrcond, ':mo', 'LineWidth',2, 'MarkerSize',8)
ylim([0,0.1])
grid on


% Add text labels beside each data point
for i = 1:length(sensor_number)
    if sensor_number(i) == 7
        text(sensor_number(i)*3, min_reciprocal_condition_number_set_minrcond(i), sprintf('%.4f', min_reciprocal_condition_number_set_minrcond(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    else
        text(sensor_number(i)*3, min_reciprocal_condition_number_set_minrcond(i), sprintf('%.4f', min_reciprocal_condition_number_set_minrcond(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end

hold on

legend('multi-objective, optimized config',  'multi-objective, larger workspace','maximize min rcond, optimized config', 'maximize min rcond, larger workspace')

sensor_number = [2;3;4;5;6;7;8;9;10];
% Set x-axis ticks to increments of 1
xticks(min(sensor_number)*3:1:max(sensor_number)*3);

xlabel("Number of axis")
ylabel("Min rcond")
title("3-axis sensors: Min rcond in different workspaces for different # of axis used")

%% 1-axis - Multiobj: Plot the minimum condition number among the 5 optimized configurations and a larger workspace
max_min_rcond_set_multiobj_1axis = [];
min_reciprocal_condition_number_set_multiobj_1axis = [];

for i = 5:18
    % Load the solution
    % Construct the file name
    fileName = sprintf('results_%d_1_axis_multiobj_2000gen_20000pop.mat', i);

    load(fileName);
    
    sens_conf = [sol];
    [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

    max_min_rcond = max(min_rcond);
    max_min_rcond_set_multiobj_1axis = [max_min_rcond_set_multiobj_1axis;max_min_rcond];

    % Get the sensor configuration that corresponds to the biggest minimum
    % condition number
    index = find(min_rcond == max(min_rcond));
    solution = sol(index,:);
    number_of_sensor = solution(end);
    solution(end) = [];
    sens_conf = reshape(solution, 7, []);
    
    % Extract position and orientation of the sensor configuration
    sens_pos = sens_conf(1:3,:);
    
    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Update the magnet configurations to a box and various orientations
    % Workspace as a box in m
    LB = [-0.05, -0.05, 0.2, deg2rad(-180), deg2rad(-90)];
    UB = [0.05, 0.05, 0.3, deg2rad(180), deg2rad(90)];
    
    % Number of samples within the box will be num_samp^3
    num_samp = 10;
    
    C = cell(1, 5);
    [C{:}] = ndgrid(linspace(0, 1, num_samp));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    magnet_conf = conf;

    % add a row of 0 for psi
    zero_row = zeros(1, size(magnet_conf,2));
    magnet_conf = [magnet_conf ; zero_row];
    
    % Get the rcond for all magnet configs in the workspace, get the
    % minimum of them
    reciprocal_condition_number_set = [];
    
    for magnet_num = 1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
        
        
        J = [];
        % Construct the Jacobian
        for sensor_num = 1:number_of_sensor
            J = [J; jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_X(sens_pos(:,sensor_num), sens_or_unitary(:,sensor_num), ...
                        magnet_pos, mu_norm, theta, phi, psi, type)];
        end
    
        sigma = svd(J);
        
        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);
    
        reciprocal_condition_number_set = [reciprocal_condition_number_set;reciprocal_condition_number];
    end
    
    min_reciprocal_condition_number_set_multiobj_1axis = [min_reciprocal_condition_number_set_multiobj_1axis; min(reciprocal_condition_number_set)];

end

for i = 1:4
    % Load the solution
    % Construct the file name
    fileName = sprintf('results_%d_1_axis_multiobj_2000gen_20000pop.mat', (i+6)*3);

    load(fileName);
    
    sens_conf = [sol];
    [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

    max_min_rcond = max(min_rcond);
    max_min_rcond_set_multiobj_1axis = [max_min_rcond_set_multiobj_1axis;max_min_rcond];

    % Get the sensor configuration that corresponds to the biggest minimum
    % condition number
    index = find(min_rcond == max(min_rcond));
    solution = sol(index,:);
    number_of_sensor = solution(end);
    solution(end) = [];
    sens_conf = reshape(solution, 7, []);
    
    % Extract position and orientation of the sensor configuration
    sens_pos = sens_conf(1:3,:);
    
    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Update the magnet configurations to a box and various orientations
    % Workspace as a box in m
    LB = [-0.05, -0.05, 0.2, deg2rad(-180), deg2rad(-90)];
    UB = [0.05, 0.05, 0.3, deg2rad(180), deg2rad(90)];
    
    % Number of samples within the box will be num_samp^3
    num_samp = 10;
    
    C = cell(1, 5);
    [C{:}] = ndgrid(linspace(0, 1, num_samp));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    magnet_conf = conf;

    % add a row of 0 for psi
    zero_row = zeros(1, size(magnet_conf,2));
    magnet_conf = [magnet_conf ; zero_row];
    
    % Get the rcond for all magnet configs in the workspace, get the
    % minimum of them
    reciprocal_condition_number_set = [];
    
    for magnet_num = 1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
        
        
        J = [];
        % Construct the Jacobian
        for sensor_num = 1:number_of_sensor
            J = [J; jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_X(sens_pos(:,sensor_num), sens_or_unitary(:,sensor_num), ...
                        magnet_pos, mu_norm, theta, phi, psi, type)];
        end
    
        sigma = svd(J);
        
        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);
    
        reciprocal_condition_number_set = [reciprocal_condition_number_set;reciprocal_condition_number];
    end
    
    min_reciprocal_condition_number_set_multiobj_1axis = [min_reciprocal_condition_number_set_multiobj_1axis; min(reciprocal_condition_number_set)];
end

figure
sensor_axis = [5;6;7;8;9;10;11;12;13;14;15;16;17;18;21;24;27;30];
plot(sensor_axis, max_min_rcond_set_multiobj_1axis, '-x', 'LineWidth',2, 'MarkerSize',8)
hold on
plot(sensor_axis, min_reciprocal_condition_number_set_multiobj_1axis, '-bx', 'LineWidth',2, 'MarkerSize',8)
ylim([0,0.1])
grid on

% Add text labels beside each data point
for i = 1:length(sensor_axis)
    text(sensor_axis(i), max_min_rcond_set_multiobj_1axis(i), sprintf('%.4f', max_min_rcond_set_multiobj_1axis(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Add text labels beside each data point
for i = 1:length(sensor_axis)
    text(sensor_axis(i), min_reciprocal_condition_number_set_multiobj_1axis(i), sprintf('%.4f', min_reciprocal_condition_number_set_multiobj_1axis(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

hold on

%% 1-axis - Min rcond: Plot the minimum condition number among the 5 optimized configurations and a larger workspace
max_min_rcond_set_minrcond_1axis = [];
min_reciprocal_condition_number_set_minrcond_1axis = [];

for i = 5:18
    % Load the solution
    % Construct the file name
    fileName = sprintf('results_%d_1_axis_minrcond_2000gen_20000pop.mat', i);

    load(fileName);
    
    sens_conf = [sol];
    [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

    max_min_rcond = max(min_rcond);
    max_min_rcond_set_minrcond_1axis = [max_min_rcond_set_minrcond_1axis;max_min_rcond];

    % Get the sensor configuration that corresponds to the biggest minimum
    % condition number
    index = find(min_rcond == max(min_rcond));
    solution = sol(index,:);
    number_of_sensor = solution(end);
    solution(end) = [];
    sens_conf = reshape(solution, 7, []);
    
    % Extract position and orientation of the sensor configuration
    sens_pos = sens_conf(1:3,:);
    
    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Update the magnet configurations to a box and various orientations
    % Workspace as a box in m
    LB = [-0.05, -0.05, 0.2, deg2rad(-180), deg2rad(-90)];
    UB = [0.05, 0.05, 0.3, deg2rad(180), deg2rad(90)];
    
    % Number of samples within the box will be num_samp^3
    num_samp = 10;
    
    C = cell(1, 5);
    [C{:}] = ndgrid(linspace(0, 1, num_samp));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    magnet_conf = conf;

    % add a row of 0 for psi
    zero_row = zeros(1, size(magnet_conf,2));
    magnet_conf = [magnet_conf ; zero_row];
    
    % Get the rcond for all magnet configs in the workspace, get the
    % minimum of them
    reciprocal_condition_number_set = [];
    
    for magnet_num = 1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
        
        
        J = [];
        % Construct the Jacobian
        for sensor_num = 1:number_of_sensor
            J = [J; jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_X(sens_pos(:,sensor_num), sens_or_unitary(:,sensor_num), ...
                        magnet_pos, mu_norm, theta, phi, psi, type)];
        end
    
        sigma = svd(J);
        
        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);
    
        reciprocal_condition_number_set = [reciprocal_condition_number_set;reciprocal_condition_number];
    end
    
    min_reciprocal_condition_number_set_minrcond_1axis = [min_reciprocal_condition_number_set_minrcond_1axis; min(reciprocal_condition_number_set)];

end

for i = 1:4
    % Load the solution
    % Construct the file name
    fileName = sprintf('results_%d_1_axis_minrcond_2000gen_20000pop.mat', (i+6)*3);

    load(fileName);
    
    sens_conf = [sol];
    [obj_rcond, min_rcond] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);

    max_min_rcond = max(min_rcond);
    max_min_rcond_set_minrcond_1axis = [max_min_rcond_set_minrcond_1axis;max_min_rcond];

    % Get the sensor configuration that corresponds to the biggest minimum
    % condition number
    index = find(min_rcond == max(min_rcond));
    solution = sol(index,:);
    number_of_sensor = solution(end);
    solution(end) = [];
    sens_conf = reshape(solution, 7, []);
    
    % Extract position and orientation of the sensor configuration
    sens_pos = sens_conf(1:3,:);
    
    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Update the magnet configurations to a box and various orientations
    % Workspace as a box in m
    LB = [-0.05, -0.05, 0.2, deg2rad(-180), deg2rad(-90)];
    UB = [0.05, 0.05, 0.3, deg2rad(180), deg2rad(90)];
    
    % Number of samples within the box will be num_samp^3
    num_samp = 10;
    
    C = cell(1, 5);
    [C{:}] = ndgrid(linspace(0, 1, num_samp));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    magnet_conf = conf;

    % add a row of 0 for psi
    zero_row = zeros(1, size(magnet_conf,2));
    magnet_conf = [magnet_conf ; zero_row];
    
    % Get the rcond for all magnet configs in the workspace, get the
    % minimum of them
    reciprocal_condition_number_set = [];
    
    for magnet_num = 1:size(magnet_conf,2)
        magnet_pos = magnet_conf(1:3,magnet_num);
        theta = magnet_conf(4, magnet_num);
        phi = magnet_conf(5, magnet_num);
        psi = magnet_conf(6, magnet_num);
        
        
        J = [];
        % Construct the Jacobian
        for sensor_num = 1:number_of_sensor
            J = [J; jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_X(sens_pos(:,sensor_num), sens_or_unitary(:,sensor_num), ...
                        magnet_pos, mu_norm, theta, phi, psi, type)];
        end
    
        sigma = svd(J);
        
        num_dof = 5;
        reciprocal_condition_number = sigma(num_dof)/sigma(1);
    
        reciprocal_condition_number_set = [reciprocal_condition_number_set;reciprocal_condition_number];
    end
    
    min_reciprocal_condition_number_set_minrcond_1axis = [min_reciprocal_condition_number_set_minrcond_1axis; min(reciprocal_condition_number_set)];
end


sensor_axis = [5;6;7;8;9;10;11;12;13;14;15;16;17;18;21;24;27;30];
plot(sensor_axis, max_min_rcond_set_minrcond_1axis, ':ro', 'LineWidth',2, 'MarkerSize',8)
hold on
plot(sensor_axis, min_reciprocal_condition_number_set_minrcond_1axis, ':mo', 'LineWidth',2, 'MarkerSize',8)
ylim([0,0.1])
grid on

% Add text labels beside each data point
for i = 1:length(sensor_axis)
    text(sensor_axis(i), max_min_rcond_set_minrcond_1axis(i), sprintf('%.4f', max_min_rcond_set_minrcond_1axis(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Add text labels beside each data point
for i = 1:length(sensor_axis)
    text(sensor_axis(i), min_reciprocal_condition_number_set_minrcond_1axis(i), sprintf('%.4f', min_reciprocal_condition_number_set_minrcond_1axis(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

hold on

xlim([5,30])
xticks(5:1:30);
xlabel("Number of axis")
ylabel("Min rcond")
title("1-axis sensors: Min rcond in different workspaces for different # of axis used")

legend('multi-objective, optimized config',  'multi-objective, larger workspace','maximize min rcond, optimized config', 'maximize min rcond, larger workspace')



%% Function
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

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame, euler angle
function J = jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_XYX(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, psi, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    Rx = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];
    
    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];

    Rx_2 = [1 0 0;                  % rotate about x-axis
    0 cos(psi) -sin(psi);
    0 sin(psi) cos(psi)];

    mu_world = Rx * Ry * Rx_2 * [0;0;mu_norm];
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
    dmu_dangle = mu_norm*[0 cos(phi)*cos(psi) -sin(phi)*sin(psi); 
                          sin(theta)*sin(psi)-cos(phi)*cos(psi)*cos(theta) sin(phi)*cos(psi)*sin(theta) -cos(theta)*cos(psi)+cos(phi)*sin(psi)*sin(theta); 
                          -cos(phi)*cos(psi)*sin(theta)-sin(psi)*cos(theta) -sin(phi)*cos(psi)*cos(theta) -cos(phi)*sin(psi)*cos(theta)-cos(psi)*sin(theta)];

    J_angle = (mu0/(4*pi))*(1/(norm(r_world)^3))*(3*(r_world_hat*r_world_hat.')-eye(3))*dmu_dangle;

    J = [J_position, J_angle];

    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
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
            psi = magnet_conf(6, magnet_num);
            
            J = [];
            for i=1:sens_num
                J = [J;jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle_XYX(sens_pos(:,i), sens_or_unitary(:,i), ...
                    magnet_pos, mu_norm, theta, phi, psi, type)];
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