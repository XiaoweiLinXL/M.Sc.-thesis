close all
clear
clc

%% Constant
scale = 1;
sph_dia = 3.175e-3*[0; 0; 1]/scale;
Volumn = (4/3)*pi*(norm(sph_dia)/2)^3;
B_r = 1.32;


%% Magnet configurations with orientation
% Workspace as a plane in unit
x = [-0.05, -0.025, 0, -0.025, 0.05]/scale;
y = [-0.05, -0.025, 0, 0.025, 0.05]/scale;
z = [0.05, 0.075, 0.10, 0.125, 0.15]/scale;
theta = linspace(-pi,pi,11);
theta(end) = [];
phi = linspace(-pi/2,pi/2,10);
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
% magnet_conf = magnet_conf(:,11046);


%% Construct the ground truth magnetic field
% Sensor axis location
P_axis = [-0.12,-0.12,0;-0.12,-0.12,0;-0.12,-0.12,0;
          -0.12,0.12,0;-0.12,0.12,0;-0.12,0.12,0;
          0.12,-0.12,0;0.12,-0.12,0;0.12,-0.12,0;
          0.12,0.12,0;0.12,0.12,0;0.12,0.12,0].';

% Sensor axis orientation
V_axis = [1,0,0; 0 1 0; 0 0 1; 
          1,0,0; 0 1 0; 0 0 1;
          1,0,0; 0 1 0; 0 0 1;
          1,0,0; 0 1 0; 0 0 1].';

B_star_axis_collection = [];
p_star_mu_star_collection = [];

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

    R_star = Rx_1*Ry*Rx_2;

    mu_hat_star = R_star * [0;0;1];

    B_star_axis = mag_field_collection(P_axis, V_axis, p_star, mu_hat_star, B_r, Volumn)*1e6;
    B_star_axis_collection = [B_star_axis_collection; B_star_axis.'];

    p_star_mu_star_collection = [p_star_mu_star_collection, [p_star;mu_hat_star]];
end
B_star_axis_collection = B_star_axis_collection.';

%%
residual_p_collection = [];
residual_R_collection = [];
residual_B_collection = [];
PM_opt_collection = [];
tic
for i = 1:size(B_star_axis_collection, 2)
    B_star_axis = B_star_axis_collection(:,i);
    P_magnet_0 = [0;0;0.10];
    Mu_hat_magnet_0 = [0;0;1];
    
    PM_0 = [P_magnet_0(:); Mu_hat_magnet_0(:)];
    
    residual_0 = residual_fun(P_axis, V_axis, ...
                              reshape(PM_0(1:3), [3,1]), reshape(PM_0(4:6), [3,1]), ...
                              Volumn, B_r, B_star_axis);
    
    
    % Define the objective function for lsqnonlin
    objective = @(PM) residual_fun(P_axis, V_axis, ...
                                   reshape(PM(1:3), [3,1]), reshape(PM(4:6), [3,1]), ...
                                   Volumn, B_r, B_star_axis);
    
    % Define the constraint function
    nonlcon = @(PM) constraints_fixed(PM);
    
    % Set optimization options
    options = optimoptions('lsqnonlin', 'Display', 'none', ...
                           'Algorithm', 'interior-point', ...
                           'OptimalityTolerance', 1e-30, ...
                           'FunctionTolerance',1e-30, ...
                           'StepTolerance', 1e-7, ...
                           'MaxFunctionEvaluations', 1e4);
    % options = optimoptions('fminimax', 'Display', 'none', ...
    %                        'MaxFunctionEvaluations',1e6, ...
    %                        'MaxIterations',1e4, ...
    %                        'FunctionTolerance', 1e-30, ...
    %                        'OptimalityTolerance', 1e-30);


    % Bounds
    lb = [-0.05, -0.05, 0.05, -inf, -inf, -inf];
    ub = [0.05, 0.05, 0.15, inf, inf, inf];
 
    
    % Call lsqnonlin
    [PM_opt,resnorm,residual,exitflag,output] = lsqnonlin(objective, PM_0, lb, ub, [], [], [], [], ...
                                                          nonlcon, options);

    % Put the solution into the list
    PM_opt_collection = [PM_opt_collection, PM_opt];
    
    % Extract optimized P_magnet and Mu_hat_magnet
    P_magnet_opt = reshape(PM_opt(1:3), [3, 1]);
    Mu_hat_magnet_opt = reshape(PM_opt(4:6), [3, 1]);

    residual_p_collection = [residual_p_collection, ...
                             norm(P_magnet_opt-p_star_mu_star_collection(1:3,i))];

    residual_R_collection = [residual_R_collection, ...
                             norm(cross(Mu_hat_magnet_opt, p_star_mu_star_collection(4:6,i)))];
    
    residual_B_collection = [residual_B_collection, ...
                             residual_fun(P_axis, V_axis, ...
                                          reshape(PM_opt(1:3), [3,1]), reshape(PM_opt(4:6), [3,1]), ...
                                          Volumn, B_r, B_star_axis)]; 
end
toc

sum(residual_p_collection>0.001)
sum(residual_R_collection>0.001)
% options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt', ...
% 'FunctionTolerance', 1e-30, 'OptimalityTolerance', 1e-30, 'StepTolerance', ...
% 1e-30, 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6);
% fun = @(x) residual_fun(x, sens_pos_collection, B_star, B_r, Volumn);
% [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[],[],options);


%% Constraint
function [c, ceq] = constraints_fixed(PM)
    % Extract Mu_hat_magnet from PM
    idx_Mu_hat_magnet = 4:6;
    Mu_hat_magnet = reshape(PM(idx_Mu_hat_magnet), [3,1]);

    % No inequality constraints (c = [])
    c = [];

    % Equality constraints for Mu_hat_magnet norms
    ceq = zeros(1, 1);
    for i = 1:1
        ceq(i) = norm(Mu_hat_magnet(:, i)) - 1;
    end
end


%%
function residual = residual_fun(P_axis, V_axis, p_magnet, mu_hat_magnet, Volumn, B_r, B_star)
%     B_star = reshape(B_star, [], 1);
    B_estimate = mag_field_collection(P_axis, V_axis, p_magnet, mu_hat_magnet, B_r, Volumn)*1e6;
%     B_collection = reshape(B_collection, [], 1) * 1e6;

    residual = abs(B_star - B_estimate);
end


function B_collection = mag_field_collection(P_axis, V_axis, p_magnet, mu_hat_magnet, B_r, Volumn)
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
        B_axis = mag_field_axis(p_axis, v_axis, p_magnet, mu_hat_magnet, B_r, Volumn);
        
        % Put that into the collection
        B_collection = [B_collection; B_axis];        
    end
end

function B_axis = mag_field_axis(p_axis, v_axis, p_magnet, mu_hat_magnet, B_r, Volumn)
    r = p_axis-p_magnet;
    r_hat = r/norm(r);
    B = (B_r*Volumn)/(4*pi*(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*mu_hat_magnet;

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