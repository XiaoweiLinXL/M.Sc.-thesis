%% Plot magnetic field
% Magnetic dipole moment
mu0 =  4*pi*1e-7; % air permeability
sph_dia = 3.175e-3*[0; 0; 1]; % sphere diameter
sph_dip = 13.2*4/3*pi*(norm(sph_dia)/2)^3/mu0; % spheres magnetic dipole
mu = sph_dip*repmat([0; 0; 1], [1, 1]); % spheres dipole
magnet_pos = [0;0;0];

% Number of points to generate along each dimension
numPointsPerDimension = 5;

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

% % Plot the points with smaller size
% scatter3(x_centers, y_centers, z_centers, 'filled', 'MarkerFaceColor', 'b', 'SizeData', 1);
% 
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('Points Equally Spaced Inside the Cube');
% axis equal;

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
B_x = B(1,:) ./ sqrt(B(1,:).^2+B(2,:).^2+B(3,:).^2)
B_y = B(2,:) ./ sqrt(B(1,:).^2+B(2,:).^2+B(3,:).^2)
B_z = B(3,:) ./ sqrt(B(1,:).^2+B(2,:).^2+B(3,:).^2)

quiver3(x, y, z, B_x, B_y, B_z);
hold on

% Draw the magnet
drawCylinder(magnet_pos,[1;0;0;0],1e-2,2e-2)
%%
% Example data
x = [1, 2, 3]; % x-coordinates
y = [2, 3, 1]; % y-coordinates
z = [3, 1, 2]; % z-coordinates
u = [0.5, -0.5, 0]; % x-components of vectors
v = [0, 0.5, -0.5]; % y-components of vectors
w = [-0.5, 0, 0.5]; % z-components of vectors

% Plot points
scatter3(x, y, z, 'filled', 'MarkerFaceColor', 'b');
hold on;

% Plot vectors
quiver3(x, y, z, u, v, w, 'Color', 'r');



xlabel('x');
ylabel('y');
zlabel('z');
title('Vectors Associated with 3D Points');
axis equal;
grid on;
hold off;

%% 
mag_field([0;1;0],[0;0;1])

%% Functions
% Magnetic field
function B = mag_field(r, mu)
% Compute magnetic field at location p for a magnet of magnetic dipole
% mu
    mu0 = 4*pi*1e-7;
    hatp = r/norm(r);

    B = mu0/(4*pi*norm(r)^3)*(3*(hatp*hatp.') - eye(3))*mu;
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






