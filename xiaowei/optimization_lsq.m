close all
clear 
clc

%% Parameters
sph_dia = 3.175e-3*[0; 0; 1]; % sphere diameter
sph_num = 10; % number of spheres
sens_noise = 0.5e-9; % RMS sensor noise 

%% Constants
Id = eye(3); % 3D identity matrix
delta = 1e-7; % derivative increment
exp_ord = 10; % order of the exponential map

% Plot
colors = [[0, 0.4470, 0.7410];
    [0.8500, 0.3250, 0.0980];
    [0.9290, 0.6940, 0.1250];
    [0.4940, 0.1840, 0.5560];
    [0.4660, 0.6740, 0.1880];
    [0.25 0.80 0.54]];

% Spheres
mu0 =  4*pi*10e-7;
sph_dip = 13.2*4/3*pi*(norm(sph_dia)/2)^3/mu0; % spheres magnetic dipole
sph_mag = sph_dip*repmat([0; 0; 1], [1, sph_num]); % spheres dipole

% Sensors
sens_no_rot = [1; 0; 0; 0];
sens_rot_x = matrixToQuaternion(Rotx(pi/2)).';
sens_rot_y = matrixToQuaternion(Roty(pi/2)).';

sens_num = 12; % number of sensors
sens_pos = 40.5e-3*[Id(:, 1), Id(:, 2), -Id(:, 1), -Id(:, 2),...
    Id(:, 1), Id(:, 2), -Id(:, 1), -Id(:, 2),...
    Id(:, 1), Id(:, 2), -Id(:, 1), -Id(:, 2)]; % sensors initial positions (w.r.t. global)
sens_or = [repmat(sens_no_rot, 1, 4), repmat(sens_rot_x, 1, 4), ...
    repmat(sens_rot_y, 1, 4)]; % sensors initial orientation - quaternion (w.r.t. global)

sens_conf = [reshape([sens_pos; sens_or], [], 1)];

%% Functions
pos = @(p0, R0, gamma) direct_kin(p0, R0, gamma, sph_dia, sph_num, exp_ord);
meas = @(p0, R0, gamma, sens_pos, sens_or) out(sens_pos, sens_or, p0, R0, gamma, sph_dia, ...
    sph_num, sph_mag, exp_ord, "1D");
meas_der = @(p0, R0, gamma, sens_pos, sens_or) num_der(sens_pos, sens_or, p0, R0, gamma, sph_dia, ...
    sph_num, sph_mag, exp_ord, delta, "1D");

%% Chain configurations
% Workspace as a box
LB = [-5e-2, -5e-2, 20e-2, -pi/2, -pi/2, -pi/2, -pi/sph_num*ones];
UB = [5e-2, 5e-2, 30e-2, pi/2, pi/2, pi/2, pi/sph_num*ones];
% Number of samples within the box
num_samp = 1;

C = cell(1, 7);
[C{:}] = ndgrid(linspace(0, 1, num_samp));
C = cellfun(@(a) a(:), C, 'Uni', 0);
combo = [C{:}];
conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
    repmat(LB, [size(combo, 1), 1]);
conf = conf.';

%% Optimization
% Limits
lb = [repmat([-20e-2, -20e-2, 0e-2, -1, -1, -1, -1], 1, 12)];
ub = [repmat([20e-2, 20e-2, 4e-2, 1, 1, 1, 1], 1, 12)];

lsq_options = optimoptions('lsqnonlin', 'Display', 'final', 'Algorithm', ...
    'levenberg-marquardt', 'FunctionTolerance', 1e-12, 'OptimalityTolerance', ...
    1e-12, 'StepTolerance', 1e-12, 'MaxFunctionEvaluations', ...
    1e6, 'MaxIterations', 1e6);


% Objective function
fun = @(x, num) min_fun(conf, x, meas, meas_der, sens_noise, num);

sol_pre = sens_conf;
for sens_num = 7:7
    [sol, resnorm, residual, exitflag, output] = lsqnonlin(@(x) fun(x, sens_num), sens_conf, lb, ub, lsq_options);
    val_sol = norm(fun(sol, sens_num));
    if sens_num > 7
        if val_sol > val
            break
        else
            sol_pre = sol;
        end
    end
    val = val_sol;
end

save('results_FL1-100_lsq')


%% Functions
function [c, ceq] = nonlin_const(conf, sens_num)
    conf = reshape(conf, 7, []);
    sens_or = conf(4:7, 1:sens_num);

    c = [];
    ceq = vecnorm(sens_or) - ones(size(vecnorm(sens_or)));
end

function obj = min_fun(conf, sens_conf, meas, meas_der, noise, sens_num)
    % Extract sensors number, position and orientation
    sens_conf = reshape(sens_conf, 7, []);
    sens_pos = sens_conf(1:3, 1:sens_num);
    sens_or = sens_conf(4:7, 1:sens_num);

    % Define objective
    obj = [];
    for i = 1:size(conf, 2)
        obs = [];
        for j = 1:sens_num
            % Signal
            signal(:, j) = meas(conf(1:3, i), Rotx(conf(4, i))*Roty(conf(5, i))*...
            Rotz(conf(6, i)), conf(7, i)*[1; 0; 0], sens_pos(:, j), quaternionToMatrix(sens_or(:, j)));
            % Observability
            obs = [obs; meas_der(conf(1:3, i), Rotx(conf(4, i))*Roty(conf(5, i))*...
            Rotz(conf(6, i)), conf(7, i)*[1; 0; 0], sens_pos(:, j), quaternionToMatrix(sens_or(:, j)))];
        end
        
        % Minimum norm amongst the sensors
        if size(signal, 1) > 1
            sig_norm = vecnorm(signal);
        else
            sig_norm = abs(signal);
        end
        [~, ind] = min(sig_norm);

        % Observability condition
        sigma = svd(obs);

        % Extend objective
        obj = [obj; 0*-sig_norm(ind)/noise; -sigma(7)/sigma(1)];
    end
end

% Rotations
function R = expr(S, ord, der)
% Exponential map on SO(3)
    if nargin < 3
        der = 0;
    end
    if nargin < 2 || isempty(ord)
        ord = 10;
    end
    R = eye(3);

    for i = 1:ord
        R = R + S^i/factorial(i + der);
    end
end

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

function q = matrixToQuaternion(R)
    % Ensure the matrix is 3x3
    assert(isequal(size(R), [3 3]), 'R must be a 3x3 matrix.');

    % Preallocate the quaternion vector
    q = zeros(1, 4);
    
    % Compute the trace of the matrix
    tr = trace(R);
    
    if tr > 0
        S = sqrt(tr + 1.0) * 2; % S=4*qw 
        q(1) = 0.25 * S;
        q(2) = (R(3,2) - R(2,3)) / S;
        q(3) = (R(1,3) - R(3,1)) / S;
        q(4) = (R(2,1) - R(1,2)) / S;
    elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        S = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3)) * 2; % S=4*qx
        q(1) = (R(3,2) - R(2,3)) / S;
        q(2) = 0.25 * S;
        q(3) = (R(1,2) + R(2,1)) / S;
        q(4) = (R(1,3) + R(3,1)) / S;
    elseif R(2,2) > R(3,3)
        S = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3)) * 2; % S=4*qy
        q(1) = (R(1,3) - R(3,1)) / S;
        q(2) = (R(1,2) + R(2,1)) / S;
        q(3) = 0.25 * S;
        q(4) = (R(2,3) + R(3,2)) / S;
    else
        S = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2)) * 2; % S=4*qz
        q(1) = (R(2,1) - R(1,2)) / S;
        q(2) = (R(1,3) + R(3,1)) / S;
        q(3) = (R(2,3) + R(3,2)) / S;
        q(4) = 0.25 * S;
    end
end


function R = Rotz(angle)
    R = [cos(angle), -sin(angle), 0;
        sin(angle), cos(angle), 0;
        0, 0, 1];
end

function R = Roty(angle)
    R = [cos(angle), 0, sin(angle);
        0, 1, 0;
        -sin(angle), 0, cos(angle)];
end

function R = Rotx(angle)
    R = [1, 0, 0;
        0, cos(angle), sin(angle);
        0, -sin(angle), cos(angle)];
end

% Skew symmetric matrix
function skew_symmetric_matrix = skew(v)
    Id = eye(3);
    skew_symmetric_matrix = [cross(v, Id(:,1)),cross(v, Id(:,2)),cross(v,Id(:,3))];
end

% Kinematics
function [p, R] = direct_kin(p0, R0, gamma, d, N, ord)
% Direct kinematics of the ball chain 
    S = expr(skew(gamma), ord); % deflection

    p = zeros(3, N); % position
    R = zeros(3, 3, N); % orientation

    % Initialize
    p(:, 1) = p0;
    R(:, :, 1) = R0;

    for i = 2:N
        R(:, :, i) = R(:, :, i - 1)*S;
        p(:, i) = p(:, i - 1) + R(:, :, i)*d;
    end
end

% Magnetics
function B = mag_field(p, mu)
% Compute magnetic field at location p for a magnet of magnetic dipole
% mu
    mu0 = 4*pi*1e-7;
    hatp = p/norm(p);

    B = mu0/(4*pi*norm(p)^3)*(3*(hatp*hatp.') - eye(3))*mu;
end

% Sensors
function y = out(ps, Rs, p0, R0, gamma, d, N, mu, ord, type)
    % Find the sensor output 
    [p, R] = direct_kin(p0, R0, gamma, d, N, ord);

    if type == "1D"
        y = zeros(1, size(ps, 2));
    else
        y = zeros(3, size(ps, 2));
    end

    for i = 1:size(ps, 2)
        for j = 1:N
            meas = Rs(:, :, i)*mag_field(ps(:, i) - p(:, j), ...
                R(:, :, j)*mu(:, j));
            if type == "1D"
                y(:, i) = y(:, i) + meas(3);
            else
                y(:, i) = y(:, i) + meas;
            end
        end
    end

    y = reshape(y, [], 1);
end

function H = num_der(ps, Rs, p0, R0, gamma, d, N, mu, ord, delta, type)
    
    meas = @(p, R, g) out(ps, Rs, p, R, g, d, N, mu, ord, type);
    q0 = [1; 0; 0; 0];

    Id = eye(3);
    Id4 = eye(4);
    
        for j = 1:3
            H(:, j) = (meas(p0 + delta*Id(:, j), R0,....
                gamma) - meas(p0, R0, gamma))/delta;
        end
        for j = 1:4
            H(:, j + 3) = (meas(p0, R0*quaternionToMatrix(q0 + delta*Id4(:, j)),....
                gamma) - meas(p0, R0, gamma))/delta;
        end
        
        H(:, 8) = (meas(p0, R0, gamma) - ...
            meas(p0, R0, gamma + delta*Id(:, 1)))/delta;
end
