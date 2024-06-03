load("results_5_one_axis_meter_rcond_4workspace.mat")

%% Constants
sens_dia = 0.5e-2;
sens_hi = 1e-2;

%% Plot magnet plate
LB = [-5e-2, -5e-2];
UB = [5e-2, 5e-2];

% Define the coordinates of the corners of the square
x = [0.05, 0.05, -0.05, -0.05, 0.05];
y = [0.05, -0.05, -0.05, 0.05, 0.05];
z = [0.2, 0.2, 0.2, 0.2, 0.2];

figure
% Plot the square
plot3(x, y, z, 'LineWidth', 2, 'MarkerSize', 10);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Sensor locations');
grid on;
hold on;

%% Plot sensors plate
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

%% Plot sensor config
optimization = "gamultiobj";
if optimization == "gamultiobj"
    sens_conf = sol;
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    sens_conf = reshape(sens_conf, 7, []);
    sens_pos = sens_conf(1:3, :);
    sens_or = sens_conf(4:7, :);
    
    for i = 1:sens_num
        drawCylinder(sens_pos(:, i), sens_or(:, i), sens_dia, sens_hi)
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
    
    % Draw the magnet
    %drawCylinder([0;0;0.2],[1;0;0;0], sens_dia, sens_hi)
elseif optimization == "paretosearch"
    sens_conf = sol;
    sens_num = 3;
    sens_conf = reshape(sens_conf, 7, []);
    sens_pos = sens_conf(1:3, :);
    sens_or = sens_conf(4:7, :);
    
    for i = 1:sens_num
        drawCylinder(sens_pos(:, i), sens_or(:, i), sens_dia, sens_hi)
    end
    
    % Draw the magnet
    %drawCylinder([0;0;0.2],[1;0;0;0], sens_dia, sens_hi)
end

%% Plot sensor config from the max rcond
sens_pos = sens_conf_max_rcond(1:3,:);
sens_or = sens_conf_max_rcond(4:7,:);
sens_num = sens_num_max_rcond;

for i = 1:sens_num
    drawCylinder(sens_pos(:, i), sens_or(:, i), sens_dia, sens_hi)
end

%% Plot sensor config from the max minB
sens_pos = sens_conf_max_minB(1:3,:);
sens_or = sens_conf_max_minB(4:7,:);
sens_num = sens_num_max_minB;

for i = 1:sens_num
    drawCylinder(sens_pos(:, i), sens_or(:, i), sens_dia, sens_hi)
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
    h.FaceColor = [1,0,0]; % Set the color of the cylinder to grey
    % axis equal;
    % xlabel('X');
    % ylabel('Y');
    % zlabel('Z');
    % title('3D Rotated and Translated Cylinder');
end