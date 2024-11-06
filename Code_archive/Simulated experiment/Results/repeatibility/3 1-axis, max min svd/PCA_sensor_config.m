clear all
close all
clc

% Load all the sensor config and form the data matrix
sens_conf_collection = [];
for i = 1:20
    filename = sprintf('data_%d.mat', i);
    load(filename);
    sens_num = sol(end);
    sol(end) = [];
    sens_conf = reshape(sol, 7, []);
    
    % Extract position and orientation
    sens_pos = sens_conf(1:3,:);

    sens_or = sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;

    sens_conf = [sens_pos; sens_or_unitary];

    % construct back the sens_conf
    one_sens_conf = [];
    for i = 1:sens_num
        one_sens_conf = [one_sens_conf, sens_conf(:,i).'];
    end

    sens_conf_collection = [sens_conf_collection; one_sens_conf];
end

% sens_conf_avg = mean(sens_conf_collection, 1);
% B = sens_conf_collection - sens_conf_avg;
% 
% [coeff,score,latent] = pca(B);
% 
% rcond1_x_coordinates = score(1:20,1);
% rcond1_y_coordinates = score(1:20,2);
% 
% % rcond0_x_coordinates = score(21:40,1);
% % rcond0_y_coordinates = score(21:40,2);
% 
% plot(rcond1_x_coordinates, rcond1_y_coordinates, 'bo');
% hold on
% 
% % plot(rcond0_x_coordinates, rcond0_y_coordinates, 'bo');
% 
% xlabel('First principle component');
% ylabel('Second principle component');
% title('Data projected on the first two principle components');
% % legend('rcond=1', 'rcond=0')

%% Random configuration

random_sensor_config_collection = [];

for i = 1:20
    one_random_sensor_config = [];
    for j = 1:3
        % Generate random vector for the first two elements from (-0.2, 0.2)
        XY = -0.2 + (0.2 + 0.2) * rand(1, 2);
        
        % Generate a zero for the third element
        Z = 0;
        
        % Generate random vector for elements (4:7) from (-1,1)
        Q = -1 + (1 + 1) * rand(1, 4);
        Q = Q/norm(Q);
        
        one_random_sensor_config = [one_random_sensor_config, XY, Z, Q];
    end
    random_sensor_config_collection = [random_sensor_config_collection; one_random_sensor_config];
end

% Total config
total_config = [sens_conf_collection;random_sensor_config_collection];
total_config_avg = mean(total_config,1);
B_total_config = total_config - total_config_avg;

[coeff,score,latent] = pca(B_total_config);

rcond1_x_coordinates = score(1:20,1);
rcond1_y_coordinates = score(1:20,2);

plot(rcond1_x_coordinates, rcond1_y_coordinates, 'ro');
hold on

random_x_coordinates = score(21:40,1);
random_y_coordinates = score(21:40,2);

plot(random_x_coordinates, random_y_coordinates, 'bo');
legend("rcond = 1", "Random config")


% random_sens_conf_avg = mean(random_sensor_config_collection, 1);
% B_random_sens_conf = random_sensor_config_collection - random_sens_conf_avg;
% 
% [coeff,score,latent] = pca(B_random_sens_conf);
% 
% rcond1_x_coordinates = score(1:20,1);
% rcond1_y_coordinates = score(1:20,2);
% 
% % rcond0_x_coordinates = score(21:40,1);
% % rcond0_y_coordinates = score(21:40,2);
% 
% plot(rcond1_x_coordinates, rcond1_y_coordinates, 'ro');
% hold on
% 
% % plot(rcond0_x_coordinates, rcond0_y_coordinates, 'bo');
% legend("rcond = 1", "Random config")


