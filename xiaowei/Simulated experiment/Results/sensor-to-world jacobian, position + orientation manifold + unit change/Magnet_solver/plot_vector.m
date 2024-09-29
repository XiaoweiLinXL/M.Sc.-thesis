% Assume v is your 7D vector, where:
v(1:3) = [0, 0, 0.1];
v(4:7) = [1, 0, 0, 0];

% Extract position
position = v(1:3);

% Extract quaternion
quat = v(4:7);

% Normalize the quaternion to avoid scaling issues
quat = quat / norm(quat);

% Define a unit vector in the direction to be rotated (e.g., along x-axis)
u = [0; 0; 1];  % This can be [1; 0; 0], [0; 1; 0], or [0; 0; 1]

% Convert quaternion to a rotation matrix
rotationMatrix = quat2rotm(quat);

% Rotate the unit vector
direction = rotationMatrix * u;

% Scale down the direction to shorten the arrow
arrow_length = 0.05;  % Set the desired length of the arrow
direction = direction * arrow_length;

% Plot the vector
quiver3(position(1), position(2), position(3), direction(1), direction(2), direction(3), 'r', 'LineWidth', 2, ...
    'AutoScale', 'off');

% Set axis limits to prevent automatic rescaling
xlim([-0.05,0.05]);  % Adjust based on your data range
ylim([-0.05,0.05]);
zlim([0,0.2]);

grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Vector Visualization');

%%
% Example stacked cell array (2x1101)
stackedCellArray = {'a', 'b', 'c', 'a';
                    'x', 'y', 'z', 'x'};  % Example with repeated columns

% Initialize a logical array to keep track of unique columns
isUnique = true(1, size(stackedCellArray, 2));

% Loop through each column
for i = 2:size(stackedCellArray, 2)
    for j = 1:i-1
        if isequal(stackedCellArray(:, i), stackedCellArray(:, j))
            isUnique(i) = false;  % Mark as not unique
            break;
        end
    end
end

% Extract the unique columns
uniqueColumns = stackedCellArray(:, isUnique);

%%
% Define the radius and center of the sphere
radius = 0.05;
center = [0, 0, 0.15];

% Generate the sphere data manually
theta = linspace(0, 2*pi, 50);  % Azimuthal angle
phi = linspace(0, pi, 50);      % Polar angle

% Create a meshgrid for the angles
[Theta, Phi] = meshgrid(theta, phi);

% Parametric equations for the sphere
X = radius * sin(Phi) .* cos(Theta) + center(1);
Y = radius * sin(Phi) .* sin(Theta) + center(2);
Z = radius * cos(Phi) + center(3);

% Plot the sphere
figure;
h = surf(X, Y, Z);

% Set the sphere to a single color and make it transparent
set(h, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Improve the appearance
axis equal;       % Equal scaling for all axes
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Transparent Blue Sphere with center at (0,0,0.15) and radius of 0.05');



