%%
sigmas_max_min_svd = [];
rconds_max_min_svd = [];
Jacobians_max_min_svd = {};
for i = 1:1:5
    filename = sprintf('results_5_one_axis_max_min_svd_%d.mat', i);
    load(filename)
    sens_conf = [sol];
    [obj_rcond, min_rcond, sigma, J] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);
    sigmas_max_min_svd = [sigmas_max_min_svd, sigma];
    rconds_max_min_svd = [rconds_max_min_svd, obj_rcond];
    Jacobians_max_min_svd{end+1} = J;

end

rcond_p_max_min_svd = [];
rcond_o_max_min_svd = [];
for i = 1:1:5
    J = Jacobians_max_min_svd{i};
    Jp = J(:,1:3);
    Jo = J(:,4:5);
    sigma_p = svd(Jp);
    sigma_o = svd(Jo);
    rcond_p = sigma_p(end)/sigma_p(1);
    rcond_o = sigma_o(end)/sigma_o(1);
    rcond_p_max_min_svd = [rcond_p_max_min_svd, rcond_p];
    rcond_o_max_min_svd = [rcond_o_max_min_svd, rcond_o];
end

%%
sigmas_max_rcond = [];
rconds_max_rcond = [];
Jacobians_max_rcond = {};
for i = 1:1:5
    filename = sprintf('results_5_one_axis_max_rcond_%d.mat', i);
    load(filename)
    sens_conf = [sol];
    [obj_rcond, min_rcond, sigma, J] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);
    sigmas_max_rcond = [sigmas_max_rcond, sigma];
    rconds_max_rcond = [rconds_max_rcond, obj_rcond];
    Jacobians_max_rcond{end+1} = J;
end

rcond_p_max_rcond = [];
rcond_o_max_rcond = [];
for i = 1:1:5
    J = Jacobians_max_rcond{i};
    Jp = J(:,1:3);
    Jo = J(:,4:5);
    sigma_p = svd(Jp);
    sigma_o = svd(Jo);
    rcond_p = sigma_p(end)/sigma_p(1);
    rcond_o = sigma_o(end)/sigma_o(1);
    rcond_p_max_rcond = [rcond_p_max_rcond, rcond_p];
    rcond_o_max_rcond = [rcond_o_max_rcond, rcond_o];
end

%%
sigmas_max_min_svd_decimeter = [];
rconds_max_min_svd_decimeter = [];
Jacobians_max_min_svd_decimeter = {};
for i = 1:1:5
    filename = sprintf('results_5_one_axis_decimeter_max_min_svd_%d.mat', i);
    load(filename)
    sens_conf = [sol];
    [obj_rcond, min_rcond, sigma, J] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type);
    sigmas_max_min_svd_decimeter = [sigmas_max_min_svd_decimeter, sigma];
    rconds_max_min_svd_decimeter = [rconds_max_min_svd_decimeter, obj_rcond];
    Jacobians_max_min_svd_decimeter{end+1} = J;

end

rcond_p_max_min_svd_decimeter = [];
rcond_o_max_min_svd_decimeter = [];
for i = 1:1:5
    J = Jacobians_max_min_svd_decimeter{i};
    Jp = J(:,1:3);
    Jo = J(:,4:5);
    sigma_p = svd(Jp);
    sigma_o = svd(Jo);
    rcond_p = sigma_p(end)/sigma_p(1);
    rcond_o = sigma_o(end)/sigma_o(1);
    rcond_p_max_min_svd_decimeter = [rcond_p_max_min_svd_decimeter, rcond_p];
    rcond_o_max_min_svd_decimeter = [rcond_o_max_min_svd_decimeter, rcond_o];
end
%%

%% 
run = 1:5;
max_singular_vals_max_svd = sigmas_max_min_svd(1,:);
min_singular_vals_max_svd = sigmas_max_min_svd(5,:);
max_singular_vals_max_rcond = sigmas_max_rcond(1,:);
min_singular_vals_max_rcond = sigmas_max_rcond(5,:);

plot(run, max_singular_vals_max_svd, 'b-', 'LineWidth',2);
hold on
plot(run, min_singular_vals_max_svd, 'b--', 'LineWidth',2)
hold on
plot(run, max_singular_vals_max_rcond, 'r-', 'LineWidth',2)
hold on
plot(run, min_singular_vals_max_rcond,  'r--', 'LineWidth',2)
hold on
xticks(run);
title('max and min singular values')
xlabel('# of run')
legend('Maximize min singular value', 'Maximize min singular value', 'Maximize rcond', 'Maximize rcond')
grid on


figure
plot(run, rconds_max_min_svd, 'blue', 'LineWidth',2);
hold on
plot(run, rconds_max_rcond, 'red', 'LineWidth',2)
hold on 
% plot(run, rcond_p_max_min_svd, 'x-', 'LineWidth', 2, 'MarkerSize', 10, 'Color', 'b')
% hold on
% plot(run, rcond_o_max_min_svd, 'x--', 'LineWidth', 2, 'MarkerSize', 10, 'Color', 'b')
% hold on
% plot(run, rcond_p_max_rcond, 'o-', 'LineWidth', 2, 'MarkerSize', 10, 'Color', 'r')
% hold on
% plot(run, rcond_o_max_rcond, 'o--', 'LineWidth', 2, 'MarkerSize', 10, 'Color', 'r')
% hold on
ylim([0,1])
xticks(run);
title('rcond')
xlabel('# of run')
legend('Maximize min singular value', 'Maximize rcond', ...
       'Maximize min singular value, rcond position', 'Maximize min singular value, rcond orientation', ...
       'Maximize rcond, rcond position', 'Maximize rcond, rcond orientation')
grid on

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
function J = jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle(sens_pos, sens_or, magnet_pos, mu_norm, theta, phi, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu0 = 4*pi*1e-7;

    Rx = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];
    
    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];

    mu_world = Rx * Ry * [0;0;mu_norm];
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
    dmu_dangle = mu_norm*[0 cos(phi); -cos(theta)*cos(phi) sin(theta)*sin(phi); -sin(theta)*cos(phi) -cos(theta)*sin(phi)];

    J_angle = (mu0/(4*pi))*(1/(norm(r_world)^3))*(3*(r_world_hat*r_world_hat.')-eye(3))*dmu_dangle;

    J = [J_position, J_angle];

    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end

%% Evaluation function
function [obj_rcond, min_rcond, sigma, J] = evaluate_with_orientation(sens_conf, magnet_conf, mu_norm, type)
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
            
            J = [];
            for i=1:sens_num
                J = [J;jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle(sens_pos(:,i), sens_or_unitary(:,i), ...
                    magnet_pos, mu_norm, theta, phi, type)];
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