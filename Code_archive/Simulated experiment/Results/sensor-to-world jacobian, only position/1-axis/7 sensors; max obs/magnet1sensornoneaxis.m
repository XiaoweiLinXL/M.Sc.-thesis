%%
close all
clear
clc
%%
% Spheres
mu0 =  4*pi*1e-7; % air permeability
sph_dia = 3.175e-3*[0; 0; 1]; % sphere diameter
sph_dip = 13.2*4/3*pi*(norm(sph_dia)/2)^3/mu0; % spheres magnetic dipole
mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole

delta = 1e-7;
type = "1D";

%% Magnet configurations 1
% Workspace as a box in m
LB = [-0.05, -0.05, 0.3];
UB = [0.05, 0.05, 0.4];

% Number of samples within the box will be num_samp^3
num_samp = 0;

C = cell(1, 3);
[C{:}] = ndgrid(linspace(0, 1, num_samp));
C = cellfun(@(a) a(:), C, 'Uni', 0);
combo = [C{:}];
conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
    repmat(LB, [size(combo, 1), 1]);
conf = conf.';

magnet_conf = conf;

%% Magnet configurations 2
% Workspace as a plane in m
LB = [-0.05, -0.05];
UB = [0.05, 0.05];

% Number of samples within the box will be num_samp^3
num_samp = 2;

C = cell(1, 2);
[C{:}] = ndgrid(linspace(0, 1, num_samp));
C = cellfun(@(a) a(:), C, 'Uni', 0);
combo = [C{:}];
conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
    repmat(LB, [size(combo, 1), 1]);
conf = conf.';

magnet_conf = conf;

z_values = 0.2 * ones(1, size(magnet_conf,2));
magnet_conf = [magnet_conf;z_values];
%% Genetic algorithm
lb = [-0.2 -0.2 0.0 -1 -1 -1 -1 ...
      -0.2 -0.2 0.0 -1 -1 -1 -1 ...
      -0.2 -0.2 0.0 -1 -1 -1 -1 ...
      -0.2 -0.2 0.0 -1 -1 -1 -1 ...
      -0.2 -0.2 0.0 -1 -1 -1 -1 ...
      -0.2 -0.2 0.0 -1 -1 -1 -1 ...
      -0.2 -0.2 0.0 -1 -1 -1 -1 ...
      7];
ub = [0.2 0.2 0.0 1 1 1 1 ...
      0.2 0.2 0.0 1 1 1 1 ...
      0.2 0.2 0.0 1 1 1 1 ...
      0.2 0.2 0.0 1 1 1 1 ...
      0.2 0.2 0.0 1 1 1 1 ...
      0.2 0.2 0.0 1 1 1 1 ...
      0.2 0.2 0.0 1 1 1 1 ...
      7];
options = optimoptions(@gamultiobj,'Display','iter', 'MaxStallGeneration',7000, 'PopulationSize', 7000, 'MaxGenerations', 100000);

fun = @(x) min_fun(x, magnet_conf, mu, type);
[sol, fval, exitflag, output] = gamultiobj(fun, length(lb), [], [], [], [], lb, ub, [], length(lb), options);

save('results_7_one_axis')
%%
% Evaluate
load('results_7_one_axis_6000gen_6000pop.mat')
sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

%%
mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

max(mean_rcond)
max(min_rcond)


%% Plot the pareto plane 2D
optimal_min_rcond = min(obj_rcond, [], 1);
optimal_min_B = min(obj_B, [], 1);
pareto_plane = [optimal_min_rcond;optimal_min_B];
figure
plot(pareto_plane(1,:), pareto_plane(2,:), 'bo', 'MarkerSize', 10);
xlabel("Min rcond in workspace")
ylabel("Min B in workspace(Tesla)")
title("Pareto plane from the solver")
grid on

%% Plot the pareto plane 3D
optimal_sensor_num = sol(:,end);
optimal_min_rcond = min(obj_rcond, [], 1);
optimal_min_B = min(obj_B, [], 1);
pareto_plane = [optimal_min_rcond;optimal_min_B;optimal_sensor_num'];
figure
plot3(pareto_plane(1,:), pareto_plane(2,:), pareto_plane(3,:), 'b.', 'MarkerSize', 10);
xlabel("Min rcond in workspace")
ylabel("Min B in workspace(Tesla)")
zlabel("# of sensor")
title("Pareto plane from the solver")
grid on

%% LSQ Nonlinear 
lb = [-0.1 -0.1 0 -0.1 -0.1 0 -0.1 -0.1 0 -0.1 -0.1 0];
ub = [0.1 0.1 0.05 0.1 0.1 0.05 0.1 0.1 0.05 0.1 0.1 0.05];

% Options
options = optimoptions('lsqnonlin', 'Display', 'iter');

fun = @(x) min_fun_lsq(x, magnet_pos, mu);

x0 = lb;
[x, resnorm,x residual, exitflag, output] = lsqnonlin(fun, x0, lb, ub, options);

%% Jacobian test
Rx = [1 0 0;                  % rotate about x-axis
0 cos(pi/3.67) -sin(pi/3.67);
0 sin(pi/3.67) cos(pi/3.67)];

Ry = [cos(-pi/4.31) 0 sin(-pi/4.31);    % rotate about y-axis
0 1 0;
-sin(-pi/4.31) 0 cos(-pi/4.31)];

sensor_output([0.0;0.0;0.0], [1;1;1;1], [0.1;0.3;0.3], Rx*Ry*[0;0;0.0176], '3D')
sensor_output_from_angle([0.0;0.0;0.0], [1;1;1;1], [0.1;0.3;0.3], 0.0176, -pi/6.1, pi/3.4, '3D')

jacobian_numerial_sensor_reading_to_magnet_conf([-0.03;0.45;0.23], [1;0;0;0], [0;0;0.3], 0.0176, pi/3.67, -pi/4.31, 1e-7, '3D')
jacobian_B_sensor_to_magnet_conf([-0.03;0.45;0.23], [1;0;0;0], [0;0;0.3], Rx*Ry*[0;0;0.0176], '3D')

%% 
B0=mag_field([0;0;0]-[0;0;0.2],mu)
%%
B1=mag_field([-0.00012515017908795;-0.00105705041190808;0]-[0;0;0.2], mu)
sensor_output([-0.00012515017908795;-0.00105705041190808;0], [0.364162472868882;-0.00182529739888303;0.000587295782991377;0.011451450689857], [0;0;0.2], mu, "1D")

B2=mag_field([2.96392419399172e-05;0.00103889782674018;0]-[0;0;0.2], mu)
sensor_output([2.96392419399172e-05;0.00103889782674018;0], [0.00208156871973884;-0.313180165010259;-0.454759542920479;-0.00208969326482125], [0;0;0.2], mu, "1D")

B3=mag_field([0.000883874772034513;-0.000589177625370821;0]-[0;0;0.2], mu)
sensor_output([0.000883874772034513;-0.000589177625370821;0], [-0.000830729717806424;0.3522942226188;-0.443797711611328;0.00288414166830842], [0;0;0.2], mu, "1D")
%%
quaternionToMatrix([0;1;-0;0])
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

% Magnetic field
function B = mag_field(r, mu)
% Compute magnetic field at location p for a magnet of magnetic dipole
% mu
    mu0 = 4*pi*1e-7;
    hatp = r/norm(r);

    B = mu0/(4*pi*norm(r)^3)*(3*(hatp*hatp.') - eye(3))*mu;
end

% Sensor output based on sensor type
function output = sensor_output(sens_pos, sens_or, magnet_pos, mu, type) % Inputs are all in world frame, sens_or is a quaternion
    field_world_frame = mag_field(sens_pos-magnet_pos, mu);

    % Convert to sensor frame
    sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(sens_or); % with respect to world frame
    field_sensor_frame = sensor_rotation_matrix.' * field_world_frame;
    
    if type == "3D"
        output = field_sensor_frame;
    elseif type == "1D"
        output = field_sensor_frame(3);
    end
end

function output = sensor_output_from_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type)
    Rx = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];

    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];

    mu_world = Rx * Ry * [0;0;mu_norm];

    B_world = mag_field(sens_pos-magnet_pos, mu_world);

    % Convert to sensor frame
    sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(sens_or); % with respect to world frame
    B_world = sensor_rotation_matrix.' * B_world;

    if type == "3D"
        output = B_world;
    elseif type == "1D"
        output = B_world(3);
    end
end

% Jacobian numerical
% From B in sensor frame to magnet position in world frame
function J = jacobian_numerial(sens_pos, sens_or, magnet_pos, mu, delta, type)
    % ps: sensor position, pm: magnet position, mu: magnetic dipole moment
    % delta: increment
    % B(ps) = f(pm, ps, mu) Compute the Jacobian with respect to pm
    % J(ps) = d f(pm, ps) / d pm
    
    senser_out = sensor_output(sens_pos, sens_or, magnet_pos, mu, type);
    
    % Numerical derivative
    delta_x = magnet_pos + [delta;0;0];
    delta_y = magnet_pos + [0;delta;0];
    delta_z = magnet_pos + [0;0;delta];
    
    
    J = [(sensor_output(sens_pos, sens_or, delta_x, mu, type)-senser_out)/delta, ... 
         (sensor_output(sens_pos, sens_or, delta_y, mu, type)-senser_out)/delta, ... 
         (sensor_output(sens_pos, sens_or, delta_z, mu, type)-senser_out)/delta];
end

% Jacobian numerical
% From B in sensor frame to magnet position in world frame
function J = jacobian_numerial_sensor_reading_to_magnet_conf(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, delta, type)
    % ps: sensor position, pm: magnet position, mu: magnetic dipole moment
    % delta: increment
    % B(ps) = f(pm, ps, mu) Compute the Jacobian with respect to pm
    % J(ps) = d f(pm, ps) / d pm
    
    senser_out = sensor_output_from_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type);
    
    % Numerical derivative
    delta_x = magnet_pos + [delta;0;0];
    delta_y = magnet_pos + [0;delta;0];
    delta_z = magnet_pos + [0;0;delta];

    delta_theta = theta + delta;
    delta_phi = phi + delta;
    
    
    J = [(sensor_output_from_angle(sens_pos, sens_or, delta_x, mu_norm, theta, phi, type)-senser_out)/delta, ... 
         (sensor_output_from_angle(sens_pos, sens_or, delta_y, mu_norm, theta, phi, type)-senser_out)/delta, ... 
         (sensor_output_from_angle(sens_pos, sens_or, delta_z, mu_norm, theta, phi, type)-senser_out)/delta, ...
         (sensor_output_from_angle(sens_pos, sens_or, magnet_pos, mu_norm, delta_theta, phi, type)-senser_out)/delta, ...
         (sensor_output_from_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, delta_phi, type)-senser_out)/delta];
end

% Jacobian analytical
% Jacobian of B in world frame to magnet position in world frame
function J = jacobian_analytical_world_to_world(sens_pos, sens_or, magnet_pos, mu, type)
    mu0 = 4*pi*1e-7;
    r = sens_pos - magnet_pos;
    r_hat = r/norm(r);

    mu_hat = mu/norm(mu);

    unit_sens_or = sens_or/norm(sens_or);

    % Convert to sensor frame
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % with respect to world frame

    J = ((3*mu0*norm(mu))/(4*pi*norm(r)^4)) * ((eye(3)-5*r_hat*r_hat.')*(r_hat.'*mu_hat) + mu_hat*r_hat.' + r_hat*mu_hat.');
    
    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame
function J = jacobian_analytical_sensor_reading_to_world(sens_pos, sens_or, magnet_pos, mu, type)
    mu0 = 4*pi*1e-7;

    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r = sens_pos - magnet_pos;
    
    r = sensor_rotation_matrix.'*r; % convert from world frame to sensor frame
    r_hat = r/norm(r);
    
    mu = sensor_rotation_matrix.'*mu; % convert from world frame to sensor frame
    mu_hat = mu/norm(mu);

    J = -((3*mu0*norm(mu))/(4*pi*norm(r)^4)) * ((eye(3)-5*r_hat*r_hat.')*(r_hat.'*mu_hat) + mu_hat*r_hat.' + r_hat*mu_hat.');
    
    J = J * sensor_rotation_matrix.';
    
    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end

% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame and
% magnet orientation with respect to world frame
function J = jacobian_B_sensor_to_magnet_conf(sens_pos, sens_or, magnet_pos, mu, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Convert r, mu from world frame to sensor frame, obtain r_hat and
    % mu_hat in sensor frame
    r = sens_pos - magnet_pos;
    
    r = sensor_rotation_matrix.'*r; % convert from world frame to sensor frame
    r_hat = r/norm(r);
    
    mu_sensor = sensor_rotation_matrix.'*mu; % convert from world frame to sensor frame
    mu_sensor_hat = mu_sensor/norm(mu_sensor);

    J_position = -((3*mu0*norm(mu_sensor))/(4*pi*norm(r)^4)) * ((eye(3)-5*(r_hat*r_hat.'))*(r_hat.'*mu_sensor_hat) + mu_sensor_hat*r_hat.' + r_hat*mu_sensor_hat.');
    
    J_position = J_position * sensor_rotation_matrix.';
    
    % If 1d sensor, select just the last row
    if type == "1D"
        J_position = J_position(3, :);
    end

    % Second part of the jacobian, from sensor reading to magnet
    % orientation

    % Extract rotation information
    mu_norm = norm(mu);
    mu_normalized = mu/mu_norm; % 1st is sin(phi), 2nd is cos(theta)cos(phi), 3rd element is sin(theta)cos(phi)
    phi = asin(mu_normalized(1));
    theta = -atan(mu_normalized(2)/mu_normalized(3));

    dmu_dangle = mu_norm*[0 cos(phi); -cos(theta)*cos(phi) sin(theta)*sin(phi); -sin(theta)*cos(phi) -cos(theta)*sin(phi)];

    J_angle = (mu0/(4*pi))*(1/(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*sensor_rotation_matrix.'*dmu_dangle;

    J = [J_position, J_angle];

    
end


%% Objective function
function obj = min_fun(sens_conf, magnet_conf, mu, type)
    
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    sens_conf = reshape(sens_conf, 7, []);
    
    % Extract position and orientation
    sens_pos = sens_conf(1:3,:);

    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    % Collect the reciprocal condition number and the min magnetic field strength for each magnet configuration
    reciprocal_condition_number_set = [];
    min_magnetic_field_set = [];
    for magnet_num=1:size(magnet_conf,2)
        magnet_pos = magnet_conf(:,magnet_num);
        
        J = [];
        for i=1:sens_num
            J = [J;jacobian_analytical_sensor_reading_to_world(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type)];
        end

        sigma = svd(J);

        num_dof = size(J,2);
        reciprocal_condition_number = sigma(num_dof)/sigma(1);
        
        % Construct the minimum sensor outputs for each sensor in the list
        sensor_outputs = [];
        for i=1:sens_num
            sensor_output_one_sensor = sensor_output(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type);
            % Put the min absolute sensor reading into the list, we care
            % about the magnitude
            sensor_outputs = [sensor_outputs, min(abs(sensor_output_one_sensor))]; 
        end

        reciprocal_condition_number_set = [reciprocal_condition_number_set; reciprocal_condition_number]; % Put the rcond for this magnet conf into the list
        min_magnetic_field_set = [min_magnetic_field_set; min(sensor_outputs)]; % Put the min sensor reading of all sensors in to the list
    end

    % Minimize the negative of the min in the list -> maximize the min in
    % the list
    obj = [];
    for i = 1:size(reciprocal_condition_number_set,1)
        obj = [obj, -min(reciprocal_condition_number_set(i))];
    end
end

%% Evaluation function
function [obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type)
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
        min_sensor_reading_one_magnet_conf = [];
        for magnet_num=1:size(magnet_conf,2)

            magnet_pos = magnet_conf(:,magnet_num);
            J = [];
            for i=1:sens_num
                J = [J;jacobian_analytical_sensor_reading_to_world(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type)];
            end
        
            sigma = svd(J);
            reciprocal_condition_number = sigma(end)/sigma(1);
            
            % Construct the minimum sensor outputs for each sensor in the list
            sensor_outputs = [];
            for i=1:sens_num
                sensor_output_one_sensor = sensor_output(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type);
                % Put the min absolute sensor reading into the list, we care
                % about the magnitude
                sensor_outputs = [sensor_outputs, min(abs(sensor_output_one_sensor))]; 
            end
            
            reciprocal_number_one_magnet_conf = [reciprocal_number_one_magnet_conf;reciprocal_condition_number];
            min_sensor_reading_one_magnet_conf = [min_sensor_reading_one_magnet_conf; min(sensor_outputs)];
        end

        obj_rcond = [obj_rcond, reciprocal_number_one_magnet_conf];
        obj_B = [obj_B, min_sensor_reading_one_magnet_conf];
    end
    min_rcond = min(obj_rcond);
    min_B = min(obj_B);
end

%% Output function
function [state,options,optchanged] = gaoutfun(options,state,flag)
    persistent h1 history r
    optchanged = false;
    switch flag
        case 'init'
            h1 = figure;
            ax = gca;
            ax.XLim = [0 21];
            ax.YLim = [0 21];
            l1 = min(state.Population(:,1));
            m1 = max(state.Population(:,1));
            l2 = min(state.Population(:,2));
            m2 = max(state.Population(:,2));
            r = rectangle(ax,'Position',[l1 l2 m1-l1 m2-l2]);
            history(:,:,1) = state.Population;
            assignin('base','gapopulationhistory',history);
        case 'iter'
            % Update the history every 10 generations.
            if rem(state.Generation,10) == 0
                ss = size(history,3);
                history(:,:,ss+1) = state.Population;
                assignin('base','gapopulationhistory',history);
            end
            % Find the best objective function, and stop if it is low.
            ibest = state.Best(end);
            ibest = find(state.Score == ibest,1,'last');
            bestx = state.Population(ibest,:);
            bestf = gaintobj(bestx);
            if bestf <= 0.1
                state.StopFlag = 'y';
                disp('Got below 0.1')
            end
            % Update the plot.
            figure(h1)
            l1 = min(state.Population(:,1));
            m1 = max(state.Population(:,1));
            l2 = min(state.Population(:,2));
            m2 = max(state.Population(:,2));
            r.Position = [l1 l2 m1-l1 m2-l2];
            pause(0.1)
            % Update the fraction of mutation and crossover after 25 generations.
            if state.Generation == 25
                options.CrossoverFraction = 0.8;
                optchanged = true;
            end
        case 'done'
            % Include the final population in the history.
            ss = size(history,3);
            history(:,:,ss+1) = state.Population;
            assignin('base','gapopulationhistory',history);
        end
end

%% LSQ Object function
function obj = min_fun_lsq(sens_conf, magnet_conf, mu)
    
    sens_conf = reshape(sens_conf, 3, []);
    sens_num = 4;
    J = [];
    for i=1:sens_num
        J = [J;jacobian_analytical(sens_conf(:,i), magnet_conf, mu)];
    end

    sigma = svd(J);
    reciprocal_condition_number = sigma(end)/sigma(1);
    
    B = [];
    for i=1:sens_num
        r = sens_conf(:,sens_num) - magnet_conf;
        B = [B, norm(mag_field(r, mu))];
    end
    
    obj = [-reciprocal_condition_number]
end