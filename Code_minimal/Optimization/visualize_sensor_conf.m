%%
close all
clear
clc
%%
load("results_5_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min_single_obj.mat")
% rcond = -fval(:,1);
% sigma_min = -fval(:,2);
% 
% sorted_sigma_min = sort(sigma_min);
% mid_point_sigma_min = sorted_sigma_min(end);
% 
% index = find(abs(sigma_min-mid_point_sigma_min)<1e-10);
% 
% fval(index,1)
% fval(index,2)
% sol = sol(index(1),:)
sol = sol;

%% Sensor config
% sol = reshape(sensor_grid_config{13}(1:2,:), 1, []);
sol = [sol, 5];

%% Constants
sens_dia = 2e-2;
sens_hi = 1e-2;

%% Plot magnet plate
% Define the radius and center
r = 0.05;
centerX = 0;
centerY = 0;

% Parametric equation for the circle
theta = linspace(0, 2*pi, 100); % Divide the circle into 100 points
x = centerX + r * cos(theta);
y = centerY + r * sin(theta);

% Plot the circle
plot(x, y, 'LineWidth', 2)
grid on 
hold on

%% Plot sensors plate
% Define Upper and Lower Bounds
ub = ub*0.1;
lb = lb*0.1;
UB = ub(1); % Upper bounds [UBx, UBy]
LB = lb(1); % Lower bounds [LBx, LBy]

UB = 0.3;
LB = -0.3;

% Define the 4 corners of the square
corners = [LB, LB, 0;
           UB, LB, 0;
           UB, UB, 0;
           LB, UB, 0];
           
x = [corners(1,1), corners(2,1), corners(3,1), corners(4,1), corners(1,1)];
y = [corners(1,2), corners(2,2), corners(3,2), corners(4,2), corners(1,2)];
z = [corners(1,3), corners(2,3), corners(3,3), corners(4,3), corners(1,3)];

% Plot each edge of the box
plot3(x,y,z, 'LineWidth', 2, 'MarkerSize', 10);
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
grid on;

%% Plot sensor config
optimization = "gamultiobj";
if optimization == "gamultiobj"
    sens_conf = sol;
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    
    sens_conf = sens_conf*0.1;
    sens_conf = reshape(sens_conf, 2, []);

    % Add the line of zero as the z coordinate to the position of sensors
    z_coordinate = zeros(1, size(sens_conf, 2));
    sens_pos = [sens_conf; z_coordinate];

    % Default orientation
    default_or = [1;0;0;0];

    % Orientation for all sensors
    sens_or = repmat(default_or, 1, sens_num);
    
    for i = 1:sens_num
        radius = sens_dia/2;
        centerX = sens_pos(1,i); % Replace with the x-coordinate of the center
        centerY = sens_pos(2,i); % Replace with the y-coordinate of the center
        
        % Parametric equation for the circle
        theta = linspace(0, 2*pi, 100); % Divide the circle into 100 points
        circleX = centerX + radius * cos(theta);
        circleY = centerY + radius * sin(theta);
        
        % Plot the circle
        plot(circleX, circleY, 'LineWidth', 2, 'Color', 'black')

%         if i == 1
%             first_sensor = sens_pos(:,1);
%             second_sensor = sens_pos(:,2);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
% 
%         if i == 2
%             first_sensor = sens_pos(:,2);
%             second_sensor = sens_pos(:,3);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
% 
%         if i == 3
%             first_sensor = sens_pos(:,3);
%             second_sensor = sens_pos(:,1);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
% 
%         if i == 4
%             first_sensor = sens_pos(:,4);
%             second_sensor = sens_pos(:,5);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
% 
%         if i == 5
%             first_sensor = sens_pos(:,5);
%             second_sensor = sens_pos(:,1);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
        
    end
    view(0,90)
    legend("Magnet workspace" ,"Sensor")
    xlim([-0.5, 0.5])
    ylim([-0.5, 0.5])
    xlabel("X (m)", 'Interpreter', 'latex')
    ylabel("Y (m)", 'Interpreter', 'latex')
    axis equal
    
% %     Draw the magnet
%     drawCylinder([0;0;0.2],[1;0;0;0], sens_dia, sens_hi)

elseif optimization == "paretosearch"
    sens_conf = sol*0.1;
    sens_num = 4;
    sens_conf = reshape(sens_conf, 2, []);

    % Add the line of zero as the z coordinate to the position of sensors
    z_coordinate = zeros(1, size(sens_conf, 2));
    sens_pos = [sens_conf; z_coordinate];

    % Default orientation
    default_or = [1;0;0;0];

    % Orientation for all sensors
    sens_or = repmat(default_or, 1, sens_num);
    
    for i = 1:sens_num
        drawCylinder(sens_pos(:, i), sens_or(:, i), sens_dia, sens_hi)
        legend("Sensor")
%         if i == 1
%             first_sensor = sens_pos(:,1);
%             second_sensor = sens_pos(:,2);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
% 
%         if i == 2
%             first_sensor = sens_pos(:,2);
%             second_sensor = sens_pos(:,3);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
% 
%         if i == 3
%             first_sensor = sens_pos(:,3);
%             second_sensor = sens_pos(:,1);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
% 
%         if i == 4
%             first_sensor = sens_pos(:,4);
%             second_sensor = sens_pos(:,5);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
% 
%         if i == 5
%             first_sensor = sens_pos(:,5);
%             second_sensor = sens_pos(:,1);
%             plot3([first_sensor(1), second_sensor(1)], [first_sensor(2), second_sensor(2)], [first_sensor(3), second_sensor(3)], 'b-');
%         end
        
    end
end

% %% Plot sensor config from the max rcond
% sens_pos = sens_conf_max_rcond(1:3,:);
% sens_or = sens_conf_max_rcond(4:7,:);
% sens_num = sens_num_max_rcond;
% 
% for i = 1:sens_num
%     drawCylinder(sens_pos(:, i), sens_or(:, i), sens_dia, sens_hi)
% end
% 
% %% Plot sensor config from the max minB
% sens_pos = sens_conf_max_minB(1:3,:);
% sens_or = sens_conf_max_minB(4:7,:);
% sens_num = sens_num_max_minB;
% 
% for i = 1:sens_num
%     drawCylinder(sens_pos(:, i), sens_or(:, i), sens_dia, sens_hi)
% end

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
    surf(X, Y, Z, "FaceColor", "flat");
    % axis equal;
    % xlabel('X');
    % ylabel('Y');
    % zlabel('Z');
    % title('3D Rotated and Translated Cylinder');
end