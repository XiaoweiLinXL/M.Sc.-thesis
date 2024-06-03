close all
clear
clc

%% Constants
% Spheres
mu0 =  4*pi*1e-7; % air permeability
sph_dia = 3.175e-3*[0; 0; 1]; % sphere diameter
sph_dip = 13.2*4/3*pi*(norm(sph_dia)/2)^3/mu0; % spheres magnetic dipole
mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole
mu_norm = norm(mu);

sens_dia = 1e-2;
sens_hi = 2e-2;

delta = 1e-7;
type = "1D";

d = 0.2; % magnet height

%% Magnet configurations
% Workspace as a box
LB = [-5e-2, -5e-2];
UB = [5e-2, 5e-2];

% Number of samples within the box
num_samp = 0;

C = cell(1, 2);
[C{:}] = ndgrid(linspace(0, 1, num_samp));
C = cellfun(@(a) a(:), C, 'Uni', 0);
combo = [C{:}];
conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
    repmat(LB, [size(combo, 1), 1]);
conf = conf.';

magnet_conf = conf;

magnet_conf = [magnet_conf; d * ones(1, size(magnet_conf, 2))];

% magnet_conf = [magnet_conf, [0.0;0.0;d], [0.0;0.05;d],[0.0;-0.05;d],[0.05;0.0;d],[-0.05;0.0;d]];

% r = 0:0.001:0.1;
% for i = 1:size(r,2)
%     magnet_conf = [magnet_conf, [r(i);0.0;0.0]];
% end
magnet_conf = [magnet_conf, [0.0;0.0;0.0;0.0;0.0]];

    

%% Construct sensor configuration
% Construct the sensor position set, where each configuration consists of 3
% sensors on a equalateral triangle, with increasing radius
sens_conf = [];
r = 0:0.1:5; % radius away from the magnet

% Radius 
d = -0.2059;
for i = 1:size(r,2)
    % Construct the vector for 3 sensors on a equalateral triangle
    %       /\        x ^  
    %      /  \         |
    %     /____\   y <-- 
    
    theta = -atan2(r(i),d);
%     % 2 sensors pointing toward the magnet
%     one_sens_conf = [r(i);0;0;cos(theta/2);0;-sin(theta/2);0;
%                      -r(i);0;0;cos(theta/2);0;sin(theta/2);0;
%                      2];

%     % 3 sensors pointing toward the magnet
%     one_sens_conf = [r;0;d(i);cos(theta/2);0;-sin(theta/2);0; 
%                      -0.5*r;sqrt(3)*r/2;d(i);cos(theta/2);sqrt(3)*sin(theta/2)/2;0.5*sin(theta/2);0; 
%                      -0.5*r;-sqrt(3)*r/2;d(i);cos(theta/2);-sqrt(3)*sin(theta/2)/2;0.5*sin(theta/2);0; 
%                      3];

%     % 3 sensors pointing upward
%     one_sens_conf = [r;0;d(i);1;0;0;0; 
%                      -0.5*r;sqrt(3)*r/2;d(i);1;0;0;0; 
%                      -0.5*r;-sqrt(3)*r/2;d(i);1;0;0;0; 
%                      3];

%     % 3 sensors pointing upward, increasing distance
%     one_sens_conf = [r;0;d(i);1;0;0;0; 
%                      -0.5*r;sqrt(3)*r/2;d(i);1;0;0;0; 
%                      -0.5*r;-sqrt(3)*r/2;d(i);1;0;0;0; 
%                      3];
 
%     % 4 sensors pointing toward the magnet
%     one_sens_conf = [r;0;d(i); cos(theta/2);0;-sin(theta/2);0;
%                      0;r;d(i); cos(theta/2); sin(theta/2);0;0; 
%                      -r;0;d(i); cos(theta/2);0;sin(theta/2);0;
%                      0;-r;d(i); cos(theta/2);-sin(theta/2);0;0;
%                      4];

%     % 4 sensors pointing upward
%     one_sens_conf = [r;0;d(i); 1;0;0;0;
%                      0;r;d(i); 1;0;0;0; 
%                      -r;0;d(i); 1;0;0;0;
%                      0;-r;d(i); 1;0;0;0;
%                      4];

%         % 3 sensors pointing upward, one randomly
%         one_sens_conf = [r(i);0;0;1;0;0;0; 
%                          -0.5*r(i);sqrt(3)*r(i)/2;0;1;0;0;0; 
%                          -0.5*r(i);-sqrt(3)*r(i)/2;0;1;0;0;0;
%                          0.73;0.29;0; 1;0;0;0;
%                          4];

    % 6 sensors pointing toward the magnet
    one_sens_conf = [sqrt(3)*r(i)/2;0.5*r(i);d; cos(theta/2);0.5*sin(theta/2);-sqrt(3)*sin(theta/2)/2;0;
                     0;r(i);d; cos(theta/2);sin(theta/2);0;0;
                     -sqrt(3)*r(i)/2;0.5*r(i);d; cos(theta/2);0.5*sin(theta/2);sqrt(3)*sin(theta/2)/2;0;
                     -sqrt(3)*r(i)/2;-0.5*r(i);d; cos(theta/2);-0.5*sin(theta/2);sqrt(3)*sin(theta/2)/2;0;
                     0;-r(i);d; cos(theta/2);-sin(theta/2);0;0;
                     sqrt(3)*r(i)/2;-0.5*r(i);d; cos(theta/2);-0.5*sin(theta/2);-sqrt(3)*sin(theta/2)/2;0;
                     6];

%     % 6 sensors pointing upward
%     one_sens_conf = [sqrt(3)*r/2;0.5*r;d(i); 1;0;0;0;
%                      0;r;d(i); 1;0;0;0;
%                      -sqrt(3)*r/2;0.5*r;d(i); 1;0;0;0;
%                      -sqrt(3)*r/2;-0.5*r;d(i); 1;0;0;0;
%                      0;-r;d(i); 1;0;0;0;
%                      sqrt(3)*r/2;-0.5*r;d(i); 1;0;0;0;
%                      6];

    sens_conf = [sens_conf, one_sens_conf];

%     % Plot sensor plate
%     UB = [-0.5 -0.5 0]; % Upper bounds [UBx, UBy, UBz]
%     LB = [0.5 0.5 0]; % Lower bounds [LBx, LBy, LBz]
%     
%     % Define the 8 corners of the box
%     corners = [LB(1), LB(2), LB(3); % Lower front left corner
%                UB(1), LB(2), LB(3); % Lower front right corner
%                UB(1), UB(2), LB(3); % Lower back right corner
%                LB(1), UB(2), LB(3); % Lower back left corner
%                LB(1), LB(2), UB(3); % Upper front left corner
%                UB(1), LB(2), UB(3); % Upper front right corner
%                UB(1), UB(2), UB(3); % Upper back right corner
%                LB(1), UB(2), UB(3)]; % Upper back left corner
%     
%     % Define lines between the corners to form the box
%     edges = [1,2; 2,3; 3,4; 4,1; % Lower square
%              5,6; 6,7; 7,8; 8,5; % Upper square
%              1,5; 2,6; 3,7; 4,8]; % Vertical lines
%     
%     % Plot each edge of the box
%     figure
%     hold on;
%     for i = 1:size(edges, 1)
%         plot3(corners(edges(i,:), 1), corners(edges(i,:), 2), corners(edges(i,:), 3), 'b-');
%     end
%     axis equal;
%     xlabel('X (m)');
%     ylabel('Y (m)');
%     zlabel('Z (m)');
%     grid on;
% 
%     % Visualize the sensors
%     sens_num = one_sens_conf(end);
%     one_sens_conf(end) = [];
%     one_sens_conf = reshape(one_sens_conf, 7, []);
%     sens_pos = one_sens_conf(1:3, :);
%     sens_or = one_sens_conf(4:7, :);
%     
%     for j = 1:sens_num
%         drawCylinder(sens_pos(:, j), sens_or(:, j), sens_dia, sens_hi)
%     end
end

% d_max_rcond = [0.0203,0.0848,0.2059];
% sens_conf = [r;0;d_max_rcond;1;0;0;0; 
%              -0.5*r;sqrt(3)*r/2;d_max_rcond;1;0;0;0; 
%              -0.5*r;-sqrt(3)*r/2;d_max_rcond;1;0;0;0; 
%              3];


%% Reciprocal condition number
rcond_all_sensor_config = [];
% Construct the reciprocal numbers from each of the configuration
for i = 1:size(sens_conf, 2)
    one_sens_conf = sens_conf(:, i);

    sens_num = one_sens_conf(end);
    one_sens_conf(end) = [];
    one_sens_conf = reshape(one_sens_conf, 7, []);
    
    % Extract position and orientation
    sens_pos = one_sens_conf(1:3,:);

    sens_or = one_sens_conf(4:7,:);
    magnitudes = vecnorm(sens_or);
    sens_or_unitary = sens_or ./ magnitudes;

    % Collect the reciprocal condition number for each magnet configuration
    rcond_all_magnet_config = [];
    for magnet_num=1:size(magnet_conf,2)
        magnet_conf = magnet_conf(:,magnet_num);
        magnet_pos = magnet_conf(1:3);
        magnet_theta = magnet_conf(4);
        magnet_phi = magnet_conf(5);
        
        J = [];
        for j=1:sens_num
            J = [J;jacobian_analytical_sensor_reading_to_magnet_conf_euler_angle(sens_pos(:,j), sens_or_unitary(:,j), magnet_pos, mu_norm, magnet_theta, magnet_phi, type)];
        end
        
        sigma = svd(J);

        num_dof = size(J,2);
        reciprocal_condition_number = sigma(num_dof)/sigma(1);

        rcond_all_magnet_config = [rcond_all_magnet_config; reciprocal_condition_number];
        
%         if abs(reciprocal_condition_number-1) < 1e-2 && (abs(sens_pos(3)-0.0374)<1e-5 || abs(sens_pos(3)-0.1416)<1e-5 || abs(sens_pos(3)-0.2059)<1e-5) 
%             disp("rcond=1")
%             figure
% %             plot_magnetic_field(mu);
%             hold on
% %             plot_gradient(mu)
%             hold on
%             
%             for j = 1:sens_num
%                 drawCylinder(sens_pos(:, j), sens_or_unitary(:, j), 0.5e-2, 2e-2)
%             end
%             drawCylinder(magnet_pos,[1;0;0;0],1e-2,2e-2)
%             
%             J_max = [];
%             for i = 1:sens_num
%                 one_sens_pos = sens_pos(:,i);
%                 one_sens_or = sens_or_unitary(:,i);
%                 field_at_sensor = mag_field(sens_pos(:,i)-magnet_pos, mu);
%                 field_at_sensor_norm = field_at_sensor./norm(field_at_sensor);
% %                 quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),field_at_sensor_norm(1),field_at_sensor_norm(2),field_at_sensor_norm(3),0.1,'LineWidth',1)
% %                 text(one_sens_pos(1)+field_at_sensor_norm(1)/10,one_sens_pos(2)+field_at_sensor_norm(2)/10,one_sens_pos(3)+field_at_sensor_norm(3)/10, 'B', 'Color', 'k', 'FontSize',10)                
%                 hold on
% 
%                 jacobian_at_sensor = -jacobian_analytical_sensor_reading_to_world(sens_pos(:,i), sens_or_unitary(:,i), magnet_pos, mu, type);
%                 J_max = [J_max;jacobian_at_sensor]
%                 row_norms = sqrt(sum(jacobian_at_sensor.^2, 2));
% 
%                 jacobian_at_sensor_norm = jacobian_at_sensor./row_norms;
%                 
%                 if type == "3D"
%                     % Plot each gradient vector
%                     quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(1),jacobian_at_sensor_norm(2),jacobian_at_sensor_norm(3),0.1,'LineWidth',1)
%                     text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of Bx_s', 'Color', 'k', 'FontSize',10)
%                     hold on
%                 
%                     quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(4),jacobian_at_sensor_norm(5),jacobian_at_sensor_norm(6),0.1,'LineWidth',1)
%                     text(one_sens_pos(1)+jacobian_at_sensor_norm(4)/10,one_sens_pos(2)+jacobian_at_sensor_norm(5)/10,one_sens_pos(3)+jacobian_at_sensor_norm(6)/10, 'Gradient of By_s', 'Color', 'k', 'FontSize',10)
%                     hold on
%                 
%                     quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(7),jacobian_at_sensor_norm(8),jacobian_at_sensor_norm(9),0.1,'LineWidth',1)
%                     text(one_sens_pos(1)+jacobian_at_sensor_norm(7)/10,one_sens_pos(2)+jacobian_at_sensor_norm(8)/10,one_sens_pos(3)+jacobian_at_sensor_norm(9)/10, 'Gradient of Bz_s', 'Color', 'k', 'FontSize',10)
%                     hold on
%                 elseif type == "1D"
%                     quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(1),jacobian_at_sensor_norm(2),jacobian_at_sensor_norm(3),0.1,'LineWidth',1)
%                     text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of B_s', 'Color', 'k', 'FontSize',10)
%                     hold on
%                 end
%                 
%             end
% 
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             title('rcond = 1')
%             axis square
%             grid on
%         end
% 
%         if abs(reciprocal_condition_number-0) < 1e-2 && (abs(sens_pos(3)-0.0)<1e-5 || abs(sens_pos(3)-0.0708)<1e-5 || abs(sens_pos(3)-0.1225)<1e-5) 
%             disp("rcond=0")
%             sens_pos
%             figure
% %             plot_magnetic_field(mu);
%             hold on
% %             plot_gradient(mu)
%             hold on
%             %Visualize the sensors
%             one_sens_conf = reshape(one_sens_conf, 7, []);
%             sens_pos = one_sens_conf(1:3, :);
%             sens_or = one_sens_conf(4:7, :);
%             
%             for j = 1:sens_num
%                 drawCylinder(sens_pos(:, j), sens_or(:, j), 0.5e-2, 2e-2)
%             end
% 
%             drawCylinder(magnet_pos,[1;0;0;0],1e-2,2e-2)
%         
%             J = [];
%             for sens_num = 1:size(sens_pos,2)
%                 one_sens_pos = sens_pos(:,sens_num);
%                 field_at_sensor = mag_field(one_sens_pos-magnet_pos, mu);
%                 field_at_sensor_norm = field_at_sensor./norm(field_at_sensor);
% %                 quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),field_at_sensor_norm(1),field_at_sensor_norm(2),field_at_sensor_norm(3),0.1,'LineWidth',1)
% %                 text(one_sens_pos(1)+field_at_sensor_norm(1)/10,one_sens_pos(2)+field_at_sensor_norm(2)/10,one_sens_pos(3)+field_at_sensor_norm(3)/10, 'B', 'Color', 'k', 'FontSize',10)
%                 hold on
%                 jacobian_at_sensor = -jacobian_analytical_sensor_reading_to_world(sens_pos(:,sens_num), sens_or(:,sens_num), magnet_pos, mu, type)
%                 jacobian_at_sensor_norm = jacobian_at_sensor./norm(jacobian_at_sensor);
%                 
%                 if type == "3D"
%                     quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(1),jacobian_at_sensor_norm(2),jacobian_at_sensor_norm(3),0.1,'LineWidth',1)
%                     text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of Bx_s', 'Color', 'k', 'FontSize',10)
%                     hold on
%                 
%                     quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(4),jacobian_at_sensor_norm(5),jacobian_at_sensor_norm(6),0.1,'LineWidth',1)
%                     text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of By_s', 'Color', 'k', 'FontSize',10)
%                     hold on
%                 
%                     quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(7),jacobian_at_sensor_norm(8),jacobian_at_sensor_norm(9),0.1,'LineWidth',1)
%                     text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of Bz_s', 'Color', 'k', 'FontSize',10)
%                     hold on
%                 elseif type == "1D"
%                     quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(1),jacobian_at_sensor_norm(2),jacobian_at_sensor_norm(3),0.1,'LineWidth',1)
%                     text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of B_s', 'Color', 'k', 'FontSize',10)
%                     hold on
%                 end
% 
%                 J = [J;jacobian_at_sensor];
%             end
%             J
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             title('rcond = 0')
%             axis equal
%         end
        
    end
    rcond_all_sensor_config = [rcond_all_sensor_config, rcond_all_magnet_config];

end

%% Plot the reciprocal condition number over radius
rcond_radius = [rcond_all_sensor_config;r];
figure
plot(r, rcond_all_sensor_config,"b-");
xlabel("distance (m)")
ylabel("rcond")
title("rcond changes with radius")
hold on

%% test 3D
J1_sens = jacobian_analytical_sensor_reading_to_world([0.1;-0.3;0.05],[0.5;0.2;0.3;0.6],[0;0;0.2],mu,'3D')
R = quaternionToMatrix([0.5;0.2;0.3;0.6]/norm([0.5;0.2;0.3;0.6]))
J1_world = R * J1_sens
J1_world.'*J1_world

%% test 1D
J1_sens = jacobian_analytical_sensor_reading_to_world([-0.0228;-0.0629;-0.3],[-0.8102;-0.3507;0.2218;-0.1409],[0;0;0],mu,'1D')
J2_sens = jacobian_analytical_sensor_reading_to_world([-0.0423;-0.0493;-0.3],[0.4027;0.8613;-0.2535;0.5280],[0;0;0],mu,'1D')
J3_sens = jacobian_analytical_sensor_reading_to_world([0.0578;0.0933;-0.3],[0.4634;-0.1312;0.5843;-0.9193],[0;0;0],mu,'1D')

J = [J1_sens;J2_sens;J3_sens]
J.'*J
svd(J)





%% Plot the B and gradient of a specific sensor config
clear all
close all
clc

load("results_one_axis.mat")

sens_conf = sol;
sens_num = sens_conf(end);
sens_conf(end) = [];
one_sens_conf = reshape(sens_conf, 7, []);
one_sens_pos = one_sens_conf(1:3, :);
one_sens_or = one_sens_conf(4:7, :);
magnitudes = vecnorm(one_sens_or);
one_sens_or_unitary = one_sens_or ./ magnitudes;

% Plot sensor space
% Define Upper and Lower Bounds
UB = ub(1:3); % Upper bounds [UBx, UBy, UBz]
LB = lb(1:3); % Lower bounds [LBx, LBy, LBz]

% Define the 8 corners of the box
corners = [LB(1), LB(2), LB(3); % Lower front left corner
           UB(1), LB(2), LB(3); % Lower front right corner
           UB(1), UB(2), LB(3); % Lower back right corner
           LB(1), UB(2), LB(3); % Lower back left corner
           LB(1), LB(2), UB(3); % Upper front left corner
           UB(1), LB(2), UB(3); % Upper front right corner
           UB(1), UB(2), UB(3); % Upper back right corner
           LB(1), UB(2), UB(3)]; % Upper back left corner

% Define lines between the corners to form the box
edges = [1,2; 2,3; 3,4; 4,1; % Lower square
         5,6; 6,7; 7,8; 8,5; % Upper square
         1,5; 2,6; 3,7; 4,8]; % Vertical lines

% Plot each edge of the box
hold on;
for i = 1:size(edges, 1)
    plot3(corners(edges(i,:), 1), corners(edges(i,:), 2), corners(edges(i,:), 3), 'b-');
end
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
grid on;

% Plot the sensors
for j = 1:sens_num
    drawCylinder(one_sens_pos(:, j), one_sens_or(:, j), 0.5e-2, 2e-2)
    hold on
end
drawCylinder([0;0;0.2], [1;0;0;0], 1e-2, 2e-2)
axis equal

J = [];
for i = 1:sens_num
    jacobian_at_sensor = -jacobian_analytical_sensor_reading_to_world(one_sens_pos(:,i), one_sens_or(:,i), [0.05;0.05;0.2], mu, type)
    J = [J;jacobian_at_sensor];

    jacobian_at_sensor_norm = jacobian_at_sensor./norm(jacobian_at_sensor);
    
    if type == '3D'
        quiver3(one_sens_pos(1,i),one_sens_pos(2,i),one_sens_pos(3,i),jacobian_at_sensor_norm(1),jacobian_at_sensor_norm(2),jacobian_at_sensor_norm(3),0.1,'LineWidth',1)
        text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of Bx', 'Color', 'k', 'FontSize',10)
        hold on
    
        quiver3(one_sens_pos(1,i),one_sens_pos(2,i),one_sens_pos(3,i),jacobian_at_sensor_norm(4),jacobian_at_sensor_norm(5),jacobian_at_sensor_norm(6),0.1,'LineWidth',1)
        text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of By', 'Color', 'k', 'FontSize',10)
        hold on
    
        quiver3(one_sens_pos(1,i),one_sens_pos(2,i),one_sens_pos(3,i),jacobian_at_sensor_norm(7),jacobian_at_sensor_norm(8),jacobian_at_sensor_norm(9),0.1,'LineWidth',1)
        text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of Bz', 'Color', 'k', 'FontSize',10)
        hold on
    elseif type == '1D'
        quiver3(one_sens_pos(1),one_sens_pos(2),one_sens_pos(3),jacobian_at_sensor_norm(1),jacobian_at_sensor_norm(2),jacobian_at_sensor_norm(3),0.1,'LineWidth',1)
        text(one_sens_pos(1)+jacobian_at_sensor_norm(1)/10,one_sens_pos(2)+jacobian_at_sensor_norm(2)/10,one_sens_pos(3)+jacobian_at_sensor_norm(3)/10, 'Gradient of Bx', 'Color', 'k', 'FontSize',10)
        hold on
    end
end


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

% Jacobian analytical
% return the jacobian of sensor reading with respect to the magnet position
% in world frame
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


%% Functions
function drawCylinder(P, Q, D, H)
    % P: Center of the cylinder [x, y, z]
    % Q: Quaternion representing orientation [w, xi, yj, zk]
    % D: Diameter of the cylinder
    % H: Height of the cylinder

    % Create a unit cylinder
    [X, Y, Z] = cylinder(1, 100);
    
    % Scale the cylinder
    X = X * D/2; % Diameter/2 for radius
    Y = Y * D/2;
    Z = Z * H; % Height

    % Convert quaternion to rotation matrix
    Q_norm = Q / norm(Q); % Normalize quaternion
    w = Q_norm(1);
    x = Q_norm(2);
    y = Q_norm(3);
    z = Q_norm(4);

    R = [1 - 2*y^2 - 2*z^2,     2*x*y - 2*z*w,       2*x*z + 2*y*w;
         2*x*y + 2*z*w,         1 - 2*x^2 - 2*z^2,   2*y*z - 2*x*w;
         2*x*z - 2*y*w,         2*y*z + 2*x*w,       1 - 2*x^2 - 2*y^2];

    % Apply rotation
    for i = 1:numel(X)
        point = [X(i); Y(i); Z(i)];
        rotated_point = R * point;
        X(i) = rotated_point(1);
        Y(i) = rotated_point(2);
        Z(i) = rotated_point(3);
    end

    % Translate the cylinder to the center
    X = X + P(1);
    Y = Y + P(2);
    Z = Z + P(3);

    % Plot the cylinder
    h = surf(X, Y, Z);
    h.FaceColor = [0.5, 0.5, 0.5]; % Set the color of the cylinder to grey
    % axis equal;
    % xlabel('X');
    % ylabel('Y');
    % zlabel('Z');
    % title('3D Rotated and Translated Cylinder');
end

function plot_magnetic_field(mu)
    magnet_pos = [0;0;0];
    
    % Number of points to generate along each dimension
    numPointsPerDimension = 14;
    
    % Generate equally spaced points inside the cube
    x = linspace(-0.2, 0.2, numPointsPerDimension);
    y = linspace(-0.2, 0.2, numPointsPerDimension);
    z = linspace(-0.0, 0.5, numPointsPerDimension);
    
    % Get coordinates of the centers of the smaller cubes
    [x_centers, y_centers, z_centers] = meshgrid(x(:), y(:), z(:));
    
    % Reshape to a column vector
    x_centers = x_centers(:);
    y_centers = y_centers(:);
    z_centers = z_centers(:);
    
    points = [x_centers.'; y_centers.'; z_centers.'];
    
    B=[];
    for i=1:size(points,2)
        point = points(:,i);
        r = point - magnet_pos;
        magnetic_field = mag_field(r, mu);
        B=[B,magnetic_field];
    end
    
    % fieldsCustom(B, points)
    x=x_centers.';
    y=y_centers.';
    z=z_centers.';
    B_x = B(1,:) ./ sqrt(B(1,:).^2+B(2,:).^2+B(3,:).^2);
    B_y = B(2,:) ./ sqrt(B(1,:).^2+B(2,:).^2+B(3,:).^2);
    B_z = B(3,:) ./ sqrt(B(1,:).^2+B(2,:).^2+B(3,:).^2);
    
    quiver3(x, y, z, B_x, B_y, B_z, 1, 'LineWidth',1);
    hold on
    
    % Draw the magnet
    drawCylinder(magnet_pos,[1;0;0;0],0.5e-2,2e-2)
end

function plot_gradient(mu)

    magnet_pos = [0;0;0];
    
    % Number of points to generate along each dimension
    numPointsPerDimension = 14;
    
    % Generate equally spaced points inside the cube
    x = linspace(-0.2, 0.2, numPointsPerDimension);
    y = linspace(-0.2, 0.2, numPointsPerDimension);
    z = linspace(-0.0, 0.5, numPointsPerDimension);
    
    % Get coordinates of the centers of the smaller cubes
    [x_centers, y_centers, z_centers] = meshgrid(x(:), y(:), z(:));
    
    % Reshape to a column vector
    x_centers = x_centers(:);
    y_centers = y_centers(:);
    z_centers = z_centers(:);
    
    points = [x_centers.'; y_centers.'; z_centers.'];
    
    J=[];
    for i=1:size(points,2)
        point = points(:,i); % sensor location
        j = -jacobian_analytical_sensor_reading_to_world(point, [1;0;0;0], magnet_pos, mu, '1D');
        J = [J, j.'];
    end
    
    % fieldsCustom(B, points)
    x=x_centers.';
    y=y_centers.';
    z=z_centers.';
    J_x = J(1,:) ./ sqrt(J(1,:).^2+J(2,:).^2+J(3,:).^2);
    J_y = J(2,:) ./ sqrt(J(1,:).^2+J(2,:).^2+J(3,:).^2);
    J_z = J(3,:) ./ sqrt(J(1,:).^2+J(2,:).^2+J(3,:).^2);
    
    quiver3(x, y, z, J_x, J_y, J_z);
    hold on
    
    % Draw the magnet
    drawCylinder(magnet_pos,[1;0;0;0],1e-2,2e-2)
end

