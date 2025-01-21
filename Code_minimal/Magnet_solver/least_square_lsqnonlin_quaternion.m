close all
clear
clc

%% Constant
scale = 1;
sph_dia = 3.175e-3*[0; 0; 1]/scale;
Volumn = 24*(4/3)*pi*(norm(sph_dia)/2)^3;
B_r = 1.32;

%% Magnet config in spherical coordinate
phi_sphere = linspace(0, 2*pi, 9);
phi_sphere(end) = [];
theta_sphere = linspace(0, pi, 7);
theta_sphere(1) = [];
theta_sphere(end) = [];
theta = linspace(-pi/2+deg2rad(30),pi/2-deg2rad(30),5);
phi = linspace(-pi/2+deg2rad(30),pi/2-deg2rad(30),5);
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
% magnet_conf = magnet_conf(:,14060);


%% Construct the ground truth magnetic field
% Sensor axis location
P_axis = [0,-0.25,0;0,-0.25,0;0,-0.25,0;
          0,0.25,0;0,0.25,0;0,0.25,0;
          0.25,0,0;0.25,0,0;0.25,0,0;
          -0.25,0,0;-0.25,0,0;-0.25,0,0].';

% Sensor axis orientation
V_axis = [1,0,0; 0 1 0; 0 0 1; 
          1,0,0; 0 1 0; 0 0 1;
          1,0,0; 0 1 0; 0 0 1;
          1,0,0; 0 1 0; 0 0 1].';

B_star_axis_collection = [];
p_star_q_star_collection = [];

for i = 1:size(magnet_conf, 2)
    one_magnet_conf = magnet_conf(:, i);
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

    mu_hat_star = R_star * [0;0;1];

    quaternion_star = rotm2quat(R_star);

    B_star_axis = mag_field_collection(P_axis, V_axis, p_star, quaternion_star, B_r, Volumn)*1e6;
    B_star_axis_collection = [B_star_axis_collection; B_star_axis.'];

    p_star_q_star_collection = [p_star_q_star_collection, [p_star;quaternion_star.']];
end
B_star_axis_collection = B_star_axis_collection.';

%%
residual_p_collection = [];
residual_R_collection = [];
residual_B_collection = [];
PM_opt_collection = [];

% Initial guess
multi_start_p = {[0;0;0.15]};
multi_start_Q = {rotm2quat(roty(315)).'};

tic
for i = 1:size(B_star_axis_collection, 2)
    B_star_axis = B_star_axis_collection(:,i);

    PM_opt_collection_multistart = [];
    resnorm_opt_collection_multistart = [];
    
    % Start from multiple initial guess
    for initial_p_index = 1:length(multi_start_p)
        for initial_Q_index = 1:length(multi_start_Q)

            P_magnet_0 = multi_start_p{initial_p_index};
            Q_magnet_0 = multi_start_Q{initial_Q_index};
    
            PM_0 = [P_magnet_0(:); Q_magnet_0(:)];
            
            residual_0 = residual_fun(P_axis, V_axis, ...
                                      reshape(PM_0(1:3), [3,1]), reshape(PM_0(4:7), [4,1]), ...
                                      Volumn, B_r, B_star_axis);
            
            
            % Define the objective function for lsqnonlin
            objective = @(PM) residual_fun(P_axis, V_axis, ...
                                           reshape(PM(1:3), [3,1]), reshape(PM(4:7), [4,1]), ...
                                           Volumn, B_r, B_star_axis);
            
            % Define the constraint function
            nonlcon = @(PM) constraints_fixed(PM);
            
            % Set optimization options
            options = optimoptions('lsqnonlin', 'Display', 'none', ...
                                   'Algorithm', 'interior-point', ...
                                   'FunctionTolerance', 1e-15, ...
                                   'OptimalityTolerance', 1e-15, ...
                                   'StepTolerance', 1e-20, ...
                                   'MaxFunctionEvaluations', 1e4, ...
                                   'MaxIterations', 1e5);
            % options = optimoptions('fminimax', 'Display', 'none', ...
            %                        'MaxFunctionEvaluations',1e6, ...
            %                        'MaxIterations',1e4, ...
            %                        'FunctionTolerance', 1e-30, ...
            %                        'OptimalityTolerance', 1e-30);
            % 
            
            % Bounds
            lb = [-0.05, -0.05, 0.10, -inf, -inf, -inf, -inf];
            ub = [0.05, 0.05, 0.20, inf, inf, inf, inf];

            % lb = [];
            % ub = [];

            % Call lsqnonlin
            [PM_opt,resnorm,residual,exitflag,output] = lsqnonlin(objective, PM_0, [], [], [], [], [], [], ...
                                                                  nonlcon, options);
        
            % Put the solution into the list
            PM_opt_collection_multistart = [PM_opt_collection_multistart, PM_opt];
            resnorm_opt_collection_multistart = [resnorm_opt_collection_multistart, resnorm]; 
        end
    end

    % From the multistart solution, pick the one that has lower error
    [smallest_residual, index] = min(resnorm_opt_collection_multistart);
    PM_opt = PM_opt_collection_multistart(:, index);
    PM_opt_collection = [PM_opt_collection, PM_opt];

    % Extract optimized P_magnet and Mu_hat_magnet
    P_magnet_opt = reshape(PM_opt(1:3), [3, 1]);
    Q_magnet_opt = reshape(PM_opt(4:7), [4, 1]);

    % Count those who are Nan
    if isnan(P_magnet_opt)
        residual_p_collection = [residual_p_collection, Inf];
        residual_R_collection = [residual_R_collection, Inf];
    else
        residual_p_collection = [residual_p_collection, ...
                                 norm(P_magnet_opt-p_star_q_star_collection(1:3,i))];

        R_magnet_opt = quaternionToMatrix(Q_magnet_opt);
        R_magnet_star = quaternionToMatrix(p_star_q_star_collection(4:7,i));
    
        residual_R_collection = [residual_R_collection, ...
                                 (abs(sum(R_magnet_opt(:,3).*R_magnet_star(:,3))-1))];
        
        residual_B_collection = [residual_B_collection, ...
                                 residual_fun(P_axis, V_axis, ...
                                              reshape(PM_opt(1:3), [3,1]), reshape(PM_opt(4:7), [4,1]), ...
                                              Volumn, B_r, B_star_axis)];
    end
end
toc

sum(residual_p_collection>0.001)
sum(residual_R_collection>0.001)
% options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt', ...
% 'FunctionTolerance', 1e-30, 'OptimalityTolerance', 1e-30, 'StepTolerance', ...
% 1e-30, 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6);
% fun = @(x) residual_fun(x, sens_pos_collection, B_star, B_r, Volumn);
% [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[],[],options);

%% Extract the configs that do not converge
index_not_converge = find(residual_R_collection>0.001);

p_star_q_star_not_converge = p_star_q_star_collection(:,index_not_converge);
PM_opt_not_converge = PM_opt_collection(:,index_not_converge);

% Plot the starting config
% Plot the vector
figure
hold on
quiver3(0,0,0.1, 0,0,0.05, 'r', 'LineWidth', 2, 'MaxHeadSize', 100);

% Set axis limits to prevent automatic rescaling
xlim([-0.5,0.5]);  % Adjust based on your data range
ylim([-0.5,0.5]);
zlim([-0.2,0.2]);

% Extract the first n configs
n = 94;
p_star_q_star_not_converge_selection = p_star_q_star_not_converge(:,1:n);
PM_opt_not_converge_6 = PM_opt_not_converge(:,1:n);
colors = jet(n); % color vector

for i = 1:size(p_star_q_star_not_converge_selection,2)
    position = p_star_q_star_not_converge_selection(1:3, i);
    quat = p_star_q_star_not_converge_selection(4:7, i);
    u = [0;0;1];

    rot = quaternionToMatrix(quat);

    direction = rot*u;

    arrow_length = 0.05;
    direction = direction*arrow_length;
    
    quiver3(position(1), position(2), position(3), direction(1), direction(2), direction(3), ...
        'LineWidth', 2, 'MaxHeadSize', 100, 'Color', colors(i,:));

    position = PM_opt_not_converge_6(1:3, i);
    quat = PM_opt_not_converge_6(4:7, i);
    u = [0;0;1];

    rot = quaternionToMatrix(quat);

    direction = rot*u;

    arrow_length = 0.05;
    direction = direction*arrow_length;
    
    quiver3(position(1), position(2), position(3), direction(1), direction(2), direction(3), ...
            'LineWidth', 2, 'MaxHeadSize', 100, 'Color', colors(i,:));
end

% Define the vertices of the cube
vertices = [
    -0.05, -0.05, 0.05;  % Vertex 1
    -0.05,  0.05, 0.05;  % Vertex 2
     0.05, -0.05, 0.05;  % Vertex 3
     0.05,  0.05, 0.05;  % Vertex 4
    -0.05, -0.05, 0.15;  % Vertex 5
    -0.05,  0.05, 0.15;  % Vertex 6
     0.05, -0.05, 0.15;  % Vertex 7
     0.05,  0.05, 0.15;  % Vertex 8
];

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

% Plot the edges of the cube
for i = 1:size(edges, 1)
    plot3([vertices(edges(i, 1), 1), vertices(edges(i, 2), 1)], ...
          [vertices(edges(i, 1), 2), vertices(edges(i, 2), 2)], ...
          [vertices(edges(i, 1), 3), vertices(edges(i, 2), 3)], 'k', 'LineWidth', 2);
end
axis equal

%% Constraint
function [c, ceq] = constraints_fixed(PM)
    % Extract Mu_hat_magnet from PM
    idx_Mu_hat_magnet = 4:7;
    Q = reshape(PM(idx_Mu_hat_magnet), [4,1]);

    % No inequality constraints (c = [])
    c = [];

    % Equality constraints for Mu_hat_magnet norms
    ceq = zeros(1, 1);
    for i = 1:1
        ceq(i) = norm(Q(:, i)) - 1;
    end
end


%%
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

function residual = residual_fun(P_axis, V_axis, p_magnet, q_magnet, Volumn, B_r, B_star)
%     B_star = reshape(B_star, [], 1);
    B_estimate = mag_field_collection(P_axis, V_axis, p_magnet, q_magnet, B_r, Volumn)*1e6;
%     B_collection = reshape(B_collection, [], 1) * 1e6;

    residual = abs(B_star - B_estimate);
end


function B_collection = mag_field_collection(P_axis, V_axis, p_magnet, q_magnet, B_r, Volumn)
    % p_magnet: magnet's position (3*1)
    % mu_hat_magnet: magnet's orientation (3*1, unit vector)
    % P_axis: sensor axis position (3*3m)
    % V_axis: sensor axis orientation (3*3m)

    B_collection = [];
    % For each axis position and orientation, construct the magnetic field
    for i = 1:size(P_axis, 2)
        p_axis = P_axis(:,i);
        v_axis = V_axis(:,i);
        
        % Magnetic field for one axis
        B_axis = mag_field_axis(p_axis, v_axis, p_magnet, q_magnet, B_r, Volumn);
        
        % Put that into the collection
        B_collection = [B_collection; B_axis];        
    end
end

function B_axis = mag_field_axis(p_axis, v_axis, p_magnet, q_magnet, B_r, Volumn)
    r = p_axis-p_magnet;
    r_hat = r/norm(r);
    B = (B_r*Volumn)/(4*pi*(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*quaternionToMatrix(q_magnet)*[0;0;1];

    % Project the field onto the sensor's axis
    B_axis = dot(B, v_axis);
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