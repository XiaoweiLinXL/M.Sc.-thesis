%%
load("results_6_3axis_10000_10000.mat")

%%
% Spheres
mu0 =  4*pi*1e-7; % air permeability
sph_dia = 3.175e-3*[0; 0; 1]; % sphere diameter
sph_dip = 13.2*4/3*pi*(norm(sph_dia)/2)^3/mu0; % spheres magnetic dipole
mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole

delta = 1e-7;
type = "3D";

%% Magnet configurations
% Workspace as a box in m
LB = [-0.05, -0.05];
UB = [0.05, 0.05];

% Number of samples within the box will be num_samp^3
num_samp = 5;

C = cell(1, 2);
[C{:}] = ndgrid(linspace(0, 1, num_samp));
C = cellfun(@(a) a(:), C, 'Uni', 0);
combo = [C{:}];
conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
    repmat(LB, [size(combo, 1), 1]);
conf = conf.';

% Row to add
height = 0.2 * ones(1, size(conf, 2)); % Create a row of 0.2s

% Add the row to the existing matrix
conf = [conf; height];

magnet_conf = conf;

%% Extract sensor information

% Store the edge of pareto plane, max rcond, max minB
max_rcond = [0];
max_rcond_minB = [];
sens_conf_max_rcond = [];

max_minB = [0];
max_minB_rcond = [];
sens_conf_max_minB = [];

% Store min rcond and its corresponding minB, check if they are the
% same as max minB and its corresponding rcond
min_rcond = [1];
min_rcond_minB = [];
% sens_conf_min_rcond = []; % Do not need to store the sensor config



magnet_conf_type = "plane";

for i = 1:size(sol,1)
    sens_conf = sol(i,:);
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    sens_conf = reshape(sens_conf, 7, []);
    
    % Extract position and orientation
    sens_pos = sens_conf(1:3,:);

    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;
    
    if magnet_conf_type == "plane"
        % If there is only one plane in magnet_conf
        % For each position of magnet, get their corresponding condition number of
        % the jacobian obtain from the sensors
        rcond_one_plane = [];
        minB_one_plane = [];
        for j = 1:size(magnet_conf,2)
            one_magnet_conf = magnet_conf(:,j);
        
            % For each sensor, get their jacobians, stack them
            J = [];
            for k = 1:sens_num
                one_sensor_jacobian = jacobian_analytical_sensor_reading_to_world(sens_pos(:,k), sens_or_unitary(:,k), one_magnet_conf, mu, type);
                J = [J; one_sensor_jacobian];
            end
    
            % Get the reciprocal condition number
            sigma = svd(J);
            reciprocal_condition_number = sigma(end)/sigma(1);
            rcond_one_plane = [rcond_one_plane, reciprocal_condition_number];
            
            % Construct the minimum sensor outputs for each sensor in the list
            sensor_outputs = [];
            for l=1:sens_num
                sensor_output_one_sensor = sensor_output(sens_pos(:,l), sens_or_unitary(:,l), one_magnet_conf, mu, type);
                sensor_outputs = [sensor_outputs, min(sensor_output_one_sensor)]; % Put the min sensor reading into the list
            end
    
            minB_one_plane = [minB_one_plane, min(abs(sensor_outputs))];  
        end
    
        % Stack the reciprocal condition number with their corresponding magnet conf
        rcond_B_with_magnet_conf = [rcond_one_plane;minB_one_plane;magnet_conf];
        
        sprintf("Min rcond: %d", min(rcond_one_plane))
        sprintf("Mean rcond: %d", mean(rcond_one_plane))
        sprintf("Min B field: %d", min(minB_one_plane))
        sprintf("Mean B field: %d", mean(minB_one_plane))
    
        % Store the data for maximum min rcond (one edge of the pareto plane)
        if min(rcond_one_plane) > min(max_rcond)
            max_rcond = rcond_one_plane;
            max_rcond_minB = minB_one_plane;
            sens_conf_max_rcond = sens_conf;
            sens_num_max_rcond = sens_num;
        end
        
        % Store the data for maximum minB (the other edge of the pareto plane)
        if min(minB_one_plane)>min(max_minB)
            max_minB = minB_one_plane;
            max_minB_rcond = rcond_one_plane;
            sens_conf_max_minB = sens_conf;
            sens_num_max_minB = sens_num;
        end
    
        if min(rcond_one_plane) < min(min_rcond)
            min_rcond = rcond_one_plane;
            min_rcond_minB = minB_one_plane;
        end
   
    
    elseif magnet_conf_type == "cube"
        % If there are multiple planes in magnet_conf
        % Store the rcond and minB for the entire magnet workspace
        rcond_all_planes = [];
        minB_all_planes = [];
    
        % For each plane in magnet conf, compute and plot their rcond and min B
        for i = 1:num_samp
            one_plane_magnet_conf = magnet_conf(:,num_samp*num_samp*(i-1)+1:num_samp*num_samp*i); % Take one plane of magnet configuration
    
            % For each position of magnet, get their corresponding condition number of
            % the jacobian obtain from the sensors
            rcond_one_plane = [];
            minB_one_plane = [];
            for j = 1:size(one_plane_magnet_conf,2)
                one_magnet_conf = one_plane_magnet_conf(:,j);
            
                % For each sensor, get their jacobians, stack them
                J = [];
                for k = 1:sens_num
                    one_sensor_jacobian = jacobian_analytical_sensor_reading_to_world(sens_pos(:,k), sens_or_unitary(:,k), one_magnet_conf, mu, type);
                    J = [J; one_sensor_jacobian];
                end
    
                % Get the reciprocal condition number
                sigma = svd(J);
                reciprocal_condition_number = sigma(end)/sigma(1);
                rcond_one_plane = [rcond_one_plane, reciprocal_condition_number];
                
                % Construct the minimum sensor outputs for each sensor in the list
                sensor_outputs = [];
                for l=1:sens_num
                    sensor_output_one_sensor = sensor_output(sens_pos(:,l), sens_or_unitary(:,l), one_magnet_conf, mu, type);
                    sensor_outputs = [sensor_outputs, min(sensor_output_one_sensor)]; % Put the min sensor reading into the list
                end
        
                minB_one_plane = [minB_one_plane, min(abs(sensor_outputs))];  
            end
    
            % Store the one plane data into the workspace data
            rcond_all_planes = [rcond_all_planes; rcond_one_plane];
            minB_all_planes = [minB_all_planes; minB_one_plane];
            
            % Stack the reciprocal condition number with their corresponding magnet conf
            rcond_B_with_magnet_conf = [rcond_one_plane;minB_one_plane;one_plane_magnet_conf];
            
            sprintf("Min rcond: %d", min(rcond_one_plane))
            sprintf("Mean rcond: %d", mean(rcond_one_plane))
            sprintf("Min B field: %d", min(minB_one_plane))
            sprintf("Mean B field: %d", mean(minB_one_plane))
           
    %         % For each plane, plot its rcond and min B field
    %         x_pos_magnet = one_plane_magnet_conf(1,:);
    %         y_pos_magnet = one_plane_magnet_conf(2,:);
    %         
    %         figure
    %         % Reshape the data into a grid
    %         [X, Y] = meshgrid(unique(x_pos_magnet), unique(y_pos_magnet));
    %         Z = griddata(x_pos_magnet, y_pos_magnet, rcond_one_plane, X, Y);
    %        
    %         % Plot the surface
    %         surf(X, Y, Z);
    %         xlabel('Magnet x position');
    %         ylabel('Magnet y position');
    %         zlabel('Reciprocal Condition Number');
    %         title(['Reciprocal Condition Number for magnet at z=',num2str(0.3+(i-1)*((0.4-0.3)/(num_samp-1)))]);
    %         
    %         zlim([0,1]);
    %         
    %         figure
    %         % Reshape the data into a grid
    %         [X, Y] = meshgrid(unique(x_pos_magnet), unique(y_pos_magnet));
    %         Z = griddata(x_pos_magnet, y_pos_magnet, minB_one_plane, X, Y);
    %         
    %         % Plot the surface
    %         surf(X, Y, Z);
    %         xlabel('Magnet x position');
    %         ylabel('Magnet y position');
    %         zlabel('Min B field');
    %         title(['Min B field for magnet at z=',num2str(0.3+(i-1)*((0.4-0.3)/(num_samp-1)))]);
        end
    
        % Store the data for maximum min rcond (one edge of the pareto plane)
        if min(rcond_all_planes) > min(max_rcond)
            max_rcond = rcond_all_planes;
            max_rcond_minB = minB_all_planes;
            sens_conf_max_rcond = sens_conf;
            sens_num_max_rcond = sens_num;
        end
        
        % Store the data for maximum minB (the other edge of the pareto plane)
        if min(minB_all_planes)>min(max_minB)
            max_minB = minB_all_planes;
            max_minB_rcond = rcond_all_planes;
            sens_conf_max_minB = sens_conf;
            sens_num_max_minB = sens_num;
        end
    
        if min(rcond_one_plane) < min(min_rcond)
            min_rcond = rcond_one_plane;
            min_rcond_minB = minB_one_plane;
        end
    end
end

%%
min(max_rcond)
mean(max_rcond)
min(max_rcond_minB)
mean(max_rcond_minB)

min(max_minB)
mean(max_minB)
min(max_minB_rcond)
mean(max_minB_rcond)

min(min_rcond)
mean(min_rcond)
min(min_rcond_minB)
mean(min_rcond_minB)

%% Plot the results with the max min rcond
figure
rcond = max_rcond;
min_B = max_rcond_minB;
x_pos_magnet = magnet_conf(1,:);
y_pos_magnet = magnet_conf(2,:);

% Reshape the data into a grid
[X, Y] = meshgrid(unique(x_pos_magnet), unique(y_pos_magnet));
Z = griddata(x_pos_magnet, y_pos_magnet, rcond, X, Y);


% Plot the surface
surf(X, Y, Z);
xlabel('Magnet x position');
ylabel('Magnet y position');
zlabel('Reciprocal Condition Number');
title('Reciprocal Condition Number');

zlim([0,1]);

figure
% Reshape the data into a grid
[X, Y] = meshgrid(unique(x_pos_magnet), unique(y_pos_magnet));
Z = griddata(x_pos_magnet, y_pos_magnet, min_B, X, Y);

% Plot the surface
surf(X, Y, Z);
xlabel('Magnet x position');
ylabel('Magnet y position');
zlabel('Min B field');
title('Min B field');

%% Plot the results with the max min B
figure
rcond = max_minB_rcond;
min_B = max_minB;
x_pos_magnet = magnet_conf(1,:);
y_pos_magnet = magnet_conf(2,:);

% Reshape the data into a grid
[X, Y] = meshgrid(unique(x_pos_magnet), unique(y_pos_magnet));
Z = griddata(x_pos_magnet, y_pos_magnet, rcond, X, Y);


% Plot the surface
surf(X, Y, Z);
xlabel('Magnet x position');
ylabel('Magnet y position');
zlabel('Reciprocal Condition Number');
title('Reciprocal Condition Number');

zlim([0,1]);

figure
% Reshape the data into a grid
[X, Y] = meshgrid(unique(x_pos_magnet), unique(y_pos_magnet));
Z = griddata(x_pos_magnet, y_pos_magnet, min_B, X, Y);

% Plot the surface
surf(X, Y, Z);
xlabel('Magnet x position');
ylabel('Magnet y position');
zlabel('Min B field');
title('Min B field');



%% Test
% For each sensor, get their jacobians, stack them
J = [];
for j = 1:sens_num
    one_sensor_jacobian = jacobian_analytical(sens_pos(:,j), sens_or(:,j), [0;0;0.2], mu, type);
    J = [J; one_sensor_jacobian];
end

%% Functions
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
