%% Distance
% Sensor positions (x, y, z)
sensor_positions = [
    1, 0, 0;
    0, 1, 0;
    -1, 0, 0;
    0, -1, 0
];

% Distance measurements (d1, d2, d3, d4)
distances = [sqrt(2)+randn(size(1))*0.1; sqrt(2)+randn(size(1))*0.1; sqrt(2)+randn(size(1))*0.01; sqrt(2)+randn(size(1))*0.1];

% Objective function
objectiveFunction = @(p) sqrt((p(1) - sensor_positions(:, 1)).^2 + ...
                              (p(2) - sensor_positions(:, 2)).^2 + ...
                              (p(3) - sensor_positions(:, 3)).^2) - distances;

% Initial guess for (x0, y0, z0)
initialGuess = [1, 1, 0.5];

% Options for lsqnonlin
options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-30);
tic
% Solve for (x0, y0, z0) using non-linear least squares
estimatedPosition = lsqnonlin(objectiveFunction, initialGuess, [], [], options);
toc
x0 = estimatedPosition(1);
y0 = estimatedPosition(2);
z0 = estimatedPosition(3);

fprintf('Estimated position: (%.4f, %.4f, %.4f)\n', x0, y0, z0);


%% Full solution
% Constants
mu0 = 4 * pi * 1e-7; % Permeability of free space (T·m/A)

% Ground truth magnet position and magnetic moment
r_m_true = [0.5, 0.5, 0.5]; % meters
m_true = [1e-5, 1e-5, 1e-5]; % A·m²

% Sensor positions
sensor_positions = [
    0, 1, 0;
    1, 0, 0;
    -1, 0, 0;
    0, -1, 0
];

% Function to calculate magnetic field from dipole
dipole_field = @(r, r_m, m) (mu0 / (4 * pi * norm(r - r_m)^3)) * ...
    (3 * dot(m, (r - r_m) / norm(r - r_m)) * (r - r_m) / norm(r - r_m) - m);

% Calculate ground truth magnetic fields at sensor positions
measured_fields = zeros(size(sensor_positions));
for i = 1:size(sensor_positions, 1)
    measured_fields(i, :) = dipole_field(sensor_positions(i, :), r_m_true, m_true);
end

% Display ground truth measured fields
disp('Ground Truth Measured Magnetic Fields (T):');
disp(measured_fields);

% Initial guess for [x_m, y_m, z_m, mx, my, mz]
initial_guess = [0.4, 0.4, 0.4, 1e-5, 1e-5, 1e-5];

% Define options for the optimizer
options = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxIterations', 1000, 'FunctionTolerance', 1e-12);

% Solve the system using nonlinear least squares
params_estimated = lsqnonlin(@(params) residuals(params, sensor_positions, measured_fields, mu0), initial_guess, [], [], options);

% Extract the estimated position and magnetic moment
estimated_position = params_estimated(1:3);
estimated_moment = params_estimated(4:6);

% Display the estimated results
disp('Estimated Magnet Position:');
disp(estimated_position);
disp('Estimated Magnetic Moment:');
disp(estimated_moment);

% Display the ground truth values for comparison
disp('Ground Truth Magnet Position:');
disp(r_m_true);
disp('Ground Truth Magnetic Moment:');
disp(m_true);

%%
B = mag_field_1([0.005;0.005;0.005],[0;0;1])
norm(B)

mu0 = 4 * pi * 1e-7;
2*(mu0/(4*pi))*(norm([0;0;1])/norm([0.005;0.005;0.005])^3)

%%
mu0 = 4 * pi * 1e-7; 
sens_pos = [0.1;0.1;0.1];
magnet_pos = [0.5;0.5;0.5];
B_r = 1.32;
sph_dia = 3.175e-3*[0; 0; 1];
Volumn = (4/3)*pi*(norm(sph_dia)/2)^3;
R = rotx(30);

B = mag_field_2(sens_pos, magnet_pos, B_r, Volumn, R)
norm(B)

r = sens_pos-magnet_pos;
mu = (B_r*Volumn/mu0) *R*[0;0;1];
B1 = mag_field_1(r,mu)

%%
function B = mag_field_1(r, mu)
    % Constants
    mu0 = 4 * pi * 1e-7; % Permeability of free space (T·m/A)

    r_hat = r/norm(r);
    B = (mu0/(4*pi*norm(r)^3)) * (3*r_hat*r_hat.'-eye(3))*mu;
end

function B = mag_field_2(sens_pos, magnet_pos, B_r, Volumn, R_star)
    r = sens_pos-magnet_pos;
    r_hat = r/norm(r);
    B = (B_r*Volumn)/(4*pi*(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*R_star*[0;0;1];
end

%%
% Residuals function for least squares optimization
function res = residuals(params, sensor_positions, measured_fields, mu0)
    x_m = params(1);
    y_m = params(2);
    z_m = params(3);
    mx = params(4);
    my = params(5);
    mz = params(6);
    r_m = [x_m, y_m, z_m];
    m = [mx, my, mz];
    num_sensors = size(sensor_positions, 1);
    res = zeros(num_sensors * 3, 1);

    for i = 1:num_sensors
        r_i = sensor_positions(i, :);
        r_vec = r_m - r_i;
        r = norm(r_vec);
        r_hat = r_vec / r;
        B_expected = (mu0 / (4 * pi * r^3)) * (3 * dot(m, r_hat) * r_hat - m);
        res((i-1)*3 + 1:i*3) = measured_fields(i, :)' - B_expected';
    end
end


%% Quaternion to rotation matrix
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

%% Skew symmetric matrix
% Convert a vector to a skew symmetric matrix
function v_hat = skew(v)
    v_hat = [0 -v(3) v(2); 
             v(3) 0 -v(1);
             -v(2) v(1) 0];
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