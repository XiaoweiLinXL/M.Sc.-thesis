clear all
close all
clc

%%
% Spheres
mu0 =  4*pi*1e-7; % air permeability
sph_dia = 3.175e-3*[0; 0; 1]; % sphere diameter
sph_dip = 13.2*4/3*pi*(norm(sph_dia)/2)^3/mu0; % spheres magnetic dipole
mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole

delta = 1e-7;
type = "1D";

%% 
load("results_3_one_axis.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

figure
scatter(3 ,max(min_rcond),"black","filled")
text(3+0.1,max(min_rcond),"1st run")
hold on

%%
load("results_4_one_axis_1.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(4 ,max(min_rcond),"black","filled")
text(4+0.1,max(min_rcond),"1st run")
hold on

%%
load("results_4_one_axis_2.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(4 ,max(min_rcond),"black","filled")
text(4+0.1,max(min_rcond),"2nd run")
hold on

%%
load("results_5_one_axis.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(5 ,max(min_rcond),"black","filled")
text(5+0.1,max(min_rcond),"1st run")
hold on

%%
load("results_6_one_axis_high_1.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(6 ,max(min_rcond),"black","filled")
text(6+0.1,max(min_rcond),"1st run")
hold on
%%
load("results_6_one_axis_high_2.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(6 ,max(min_rcond),"black","filled")
text(6+0.1,max(min_rcond),"2nd run")
hold on
%%
load("results_6_one_axis_high_3.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(6 ,max(min_rcond),"black","filled")
text(6+0.1,max(min_rcond),"3rd run")
hold on
%%
load("results_6_one_axis_high_4.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(6 ,max(min_rcond),"black","filled")
text(6+0.1,max(min_rcond),"4th run")
hold on
%%
load("results_6_one_axis_high_5.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(6 ,max(min_rcond),"black","filled")
text(6+0.1,max(min_rcond),"5th run")
hold on
%%
load("results_6_one_axis_high_6.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(6 ,max(min_rcond),"black","filled")
text(6+0.1,max(min_rcond),"6th run")
hold on
%%
load("results_7_one_axis.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(7 ,max(min_rcond),"black","filled")
text(7+0.1,max(min_rcond),"1st run")
hold on
%%
load("results_8_one_axis.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(8 ,max(min_rcond),"black","filled")
text(8+0.1,max(min_rcond),"1st run")
hold on
%%
load("results_9_one_axis_high_1.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(9 ,max(min_rcond),"black","filled")
text(9+0.1,max(min_rcond),"1st run")
hold on
%%
load("results_9_one_axis_high_2.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(9 ,max(min_rcond),"black","filled")
text(9+0.1,max(min_rcond),"2nd run")
hold on
%%
load("results_9_one_axis_high_3.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(9 ,max(min_rcond),"black","filled")
text(9+0.1,max(min_rcond),"3rd run")
hold on
%%
load("results_9_one_axis_high_4.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(9 ,max(min_rcond),"black","filled")
text(9+0.1,max(min_rcond),"4th run")
hold on

%%
load("results_9_one_axis_high_5.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(9 ,max(min_rcond),"black","filled")
text(9+0.1,max(min_rcond),"5th run")
hold on

%%
load("results_9_one_axis_high_6.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(9 ,max(min_rcond),"black","filled")
text(9+0.1,max(min_rcond),"6th run")
hold on



%%
type = "3D";

load("results_1_three_axis.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(3 ,max(min_rcond),"blue","filled")
text(3+0.1, max(min_rcond), "2 3-axis sensor")
hold on



%%
load("results_2_three_axis.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(6 ,max(min_rcond),"blue","filled")
text(6+0.1, max(min_rcond), "2 3-axis sensor")
hold on

%%
load("results_3_three_axis.mat")

sens_conf = [sol];
[obj_rcond, obj_B, min_rcond, min_B] = evaluate(sens_conf, magnet_conf, mu, type);
sol
obj_rcond
obj_B
min_rcond
min_B

mean_rcond = [];
for i = 1:size(obj_rcond,2)
    mean_rcond = [mean_rcond, mean(obj_rcond(:,i))];
end

min_rcond = [];
for i = 1:size(obj_rcond, 2)
    min_rcond = [min_rcond, min(obj_rcond(:,i))];
end

scatter(9 ,max(min_rcond),"blue","filled")
text(9+0.1, max(min_rcond), "3 3-axis sensor")
hold on


xlim([0,10])
ylim([0,1])
grid on
xlabel("Number of sensors")
ylabel("Minimum rcond in the workspace")
title("Optimization results")

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