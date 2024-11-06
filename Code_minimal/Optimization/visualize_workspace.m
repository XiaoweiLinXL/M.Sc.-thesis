% phi = linspace(0, 2*pi, 9);
% phi(end) = [];
% theta = linspace(0, pi, 7);
% theta(1) = [];
% 
% [PHI, THETA] = ndgrid(phi, theta);
% 
% PHI = reshape(PHI, [], 1);
% THETA = reshape(THETA, [], 1);
% 
% ball = [PHI, THETA].';
% 
% ball = [ball, [0;0],[0;pi]];
% 
% points = [];
% 
% r = 0.05;
% for i = 1:size(ball, 2)
%     phi = ball(1, i);
%     theta = ball(2, i);
%     x = r*sin(theta)*cos(phi);
%     y = r*sin(theta)*sin(phi);
%     z = r*cos(theta);
% 
%     points = [points, [x;y;z]];
% end
% 
% x = points(1,:);
% y = points(2,:);
% z = points(3,:) + 0.15;
% 
% figure;
% hold on
% 
% plot3(x, y, z, 'o'); % 'o' specifies circle markers
% 
% grid on; % Adds grid lines for better visualization
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% axis equal
% title('3D Scatter Plot of Points');

% Create a sphere
[x, y, z] = sphere(50);  % 50 is the resolution of the sphere

% Scale the sphere to the desired radius and shift its center
radius = 0.05;
x = radius * x;
y = radius * y;
z = radius * z + 0.15;

% Plot the sphere
figure
hold on
grid on
surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 0.3 makes the sphere 30% opaque

% Set axis properties
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
% % title('Transparent Sphere with Center at (0, 0.1) and Radius 0.05');
% 
xlim([-0.075,0.075])
ylim([-0.075 0.075])
zlim([0 0.25])

% quiver3(0, 0, 0.2, 0, 0, 0.02, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% quiver3(0, 0, 0.15, 0, 0, 0.02, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% quiver3(0, 0, 0.1, 0, 0, 0.02, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% quiver3(0, 0.05, 0.15, 0, 0, 0.02, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% quiver3(0.05, 0, 0.15, 0, 0, 0.02, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% quiver3(0, -0.05, 0.15, 0, 0, 0.02, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% quiver3(-0.05, 0, 0.15, 0, 0, 0.02, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% 
% 
% 
% 
drawCylinder([0,0.05,0], [1;0;0;0], 0.5e-2, 1e-2)
drawCylinder([0.05,0,0], [1;0;0;0], 0.5e-2, 1e-2)
drawCylinder([0,-0.05,0], [1;0;0;0], 0.5e-2, 1e-2)
drawCylinder([-0.05,0,0], [1;0;0;0], 0.5e-2, 1e-2)


%%
figure
hold on
grid on;
axis equal

% Create a sphere
[x, y, z] = sphere(50);  % 50 is the resolution of the sphere

% Scale the sphere to the desired radius and shift its center
radius = 0.05;
x = radius * x;
y = radius * y;
z = radius * z + 0.15;

% % Plot the sphere
% surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 0.3 makes the sphere 30% opaque

% theta = linspace(0, 2*pi, 9);
% theta(end) = [];
% phi = linspace(0, pi, 5);
% phi(end) = [];
% psi = [0];

theta = linspace(-pi/2+deg2rad(15),pi/2-deg2rad(15),5);
phi = linspace(-pi/2+deg2rad(15),pi/2-deg2rad(15),5);


% theta = linspace(0, 2*pi, 9);
% theta(end) = [];
% phi = linspace(0, pi, 5);
% phi(end) = [];
psi = [0];

[Theta, Phi, Psi] = ndgrid(theta, phi, psi);

Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);
Psi = reshape(Psi, [], 1);

Angles = [Theta, Phi, Psi].';


v = [0;0;1];

for i = 1:size(Angles, 2)
    theta = Angles(1,i);
    phi = Angles(2,i);
    psi = Angles(3,i);

    Rx_1 = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];
    
    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];
    
    Rx_2 = [1 0 0;                  % rotate about x-axis
    0 cos(psi) -sin(psi);
    0 sin(psi) cos(psi)];

    % R_star = Rx_1 * Ry * Rx_2;
    R_star = Ry * Rx_1 * Rx_2;

    v_rotated = R_star * v *0.02;

    if i == (1+size(Angles, 2))/2
        quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'r');
    else
        quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'b');

    end


    % quiver3(0, 0, 0.10, v_rotated(1), v_rotated(2), v_rotated(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % quiver3(0, 0, 0.15, v_rotated(1), v_rotated(2), v_rotated(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % quiver3(0, 0, 0.20, v_rotated(1), v_rotated(2), v_rotated(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % quiver3(0.05, 0, 0.15, v_rotated(1), v_rotated(2), v_rotated(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % quiver3(0, 0.05, 0.15, v_rotated(1), v_rotated(2), v_rotated(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % quiver3(-0.05, 0, 0.15, v_rotated(1), v_rotated(2), v_rotated(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % quiver3(0, -0.05, 0.15, v_rotated(1), v_rotated(2), v_rotated(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % 

end
xlabel("X" ,'interpreter' , 'latex')
ylabel("Y" ,'interpreter' , 'latex')
zlabel("Z" ,'interpreter' , 'latex')


% zlim([0, 0.25])
% xlim([-0.1 0.1])
% ylim([-0.1 0.1])

% drawCylinder([0,0.10,0], [1;0;0;0], 0.5e-2, 1e-2)
% drawCylinder([0.10,0,0], [1;0;0;0], 0.5e-2, 1e-2)
% drawCylinder([0,-0.10,0], [1;0;0;0], 0.5e-2, 1e-2)
% drawCylinder([-0.10,0,0], [1;0;0;0], 0.5e-2, 1e-2)
% set(gca, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);

%% 24-side polygon
% Number of sides
n = 24;

% Angles for the vertices (0 to 2*pi)
theta = linspace(0, 2*pi, n+1);

% Define the radius (for simplicity, let's take it as 1)
r = 1;

% X and Y coordinates of the polygon vertices
x = r * cos(theta);
y = r * sin(theta);

% Plot the polygon
figure;
plot(x, y, 'b-', 'LineWidth', 2); % Blue line for the polygon
axis equal; % Equal scaling on both axes
grid on;
title('24-sided Polygon');


%% Function
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


