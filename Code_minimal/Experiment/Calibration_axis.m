close all
clear 
clc

%% Constants
cylinder_diameter = 6.35e-3; % 6.35mm, diameter of the cylindrical magnet
cylinder_length = 2*6.35e-3; % 2*6.35mm, length of the cylindrical magnet
cylinder_volume = cylinder_length*pi*(cylinder_diameter/2)^2; % Volumn of a sphere 

B_r = 1.32; % Magnetic flux density (T)

%%
% Folder path
folder = fullfile(pwd, "data_5cm");

% .mat files
matFiles = dir(fullfile(folder, '*.mat'));

B_star = [];
magnet_conf = [];

% Go through the data
for i = 1:length(matFiles)
    % Get the full path of the .mat file
    matFilePath = fullfile(folder, matFiles(i).name);

    % Load the .mat file
    data = load(matFilePath);
    
    % Process the loaded data (replace this with your actual processing code)
    disp(['Processing file: ', matFiles(i).name]);

    % Extract the filename without the extension, for the magnet
    % configuration
    [~, fileName, ~] = fileparts(matFiles(i).name);
    
    % Split the filename by the underscore "_"
    parts = strsplit(fileName, '_');
    
    % Extract the first part and the second part (if they exist)
    if length(parts) >= 2
        firstPart = parts{1};  % Part before the first underscore
        secondPart = parts{2}; % Part between the first and second underscores
    elseif length(parts) == 1
        firstPart = parts{1};  % If no underscore, entire name is the first part
        secondPart = '';       % No second part available
    else
        firstPart = '';
        secondPart = '';
    end
    
    % Display the extracted parts
    disp(['Position: ', firstPart]);
    disp(['Y orientation: ', secondPart]);

    if firstPart == "0010"
        p_m = [0;0;0.1];
    elseif  firstPart == "0015"
        p_m = [0;0;0.15];
    elseif firstPart == "0020"
        p_m = [0;0;0.2];
    elseif firstPart == "0n515"
        p_m = [0;-0.05;0.15];
    elseif firstPart == "n5015"
        p_m = [-0.05;0;0.15];
    elseif firstPart == "5015"
        p_m = [0.05;0;0.15];
    elseif firstPart == "0515"
        p_m = [0;0.05;0.15];
    end

    if secondPart == "y0"
        R_y = roty(0);
    elseif secondPart == "y45"
        R_y = roty(45);
    elseif secondPart == "y90"
        R_y = roty(90);
    elseif secondPart == "y135"
        R_y = roty(135);
    end

    
    % Go through each saved field
    for j = 1:length(data.save_element)
        B_star = [B_star, mean(data.user_data.Field_filt(:, data.save_element(j)-5: data.save_element(j)+5), 2)];

        if j == 1
            R_x = rotx(0);
        elseif j == 2
            R_x = rotx(45);
        elseif j == 3
            R_x = rotx(90);
        elseif j == 4
            R_x = rotx(135);
        elseif j == 5
            R_x = rotx(180);
        elseif j == 6
            R_x = rotx(225);
        elseif j == 7
            R_x = rotx(270);
        elseif j == 8
            R_x = rotx(315);
        end

        R = R_y * R_x;

        mu_hat = R(:, 3);
        
        magnet_conf = [magnet_conf, [p_m;mu_hat]];

    end
end

%% Calibration - Initial guess
% Initial guess for axis position
P_axis_0 = [0.05, 0, 0; 
            0.05, 0, 0; 
            0.05, 0, 0; 
            0, 0.05, 0; 
            0, 0.05, 0;
            0, 0.05, 0;
            -0.05, 0, 0; 
            -0.05, 0, 0; 
            -0.05, 0, 0; 
            0, -0.05, 0; 
            0, -0.05, 0; 
            0, -0.05, 0].';

% Initial guess for axis orientation
V_axis_0 = [1, 0, 0;
            0, 1, 0;
            0, 0, 1;
            1, 0, 0;
            0, 1, 0;
            0, 0, 1;
            1, 0, 0;
            0, 1, 0;
            0, 0, 1;
            1, 0, 0;
            0, 1, 0;
            0, 0, 1].';

PV_axis_0 = [P_axis_0(:);V_axis_0(:)]; % initial guess

num_sensor = 4; % number of sensor
num_axis = 3 * num_sensor; % number of axis

residual_before = residual_fun(P_axis_0, V_axis_0, magnet_conf(1:3,:), magnet_conf(4:6,:), cylinder_volume, B_r, B_star);

%% fminimax
% Objective function is the residual vector
objective = @(PV_axis) residual_fun(reshape(PV_axis(1:3*num_axis), [3, num_axis]), reshape(PV_axis(3*num_axis+1:end), [3, num_axis]), magnet_conf(1:3,:), magnet_conf(4:6,:), cylinder_volume, B_r, B_star);

% Define the constraint function
nonlcon = @(PV_axis) constraints(PV_axis, num_axis);

% Set optimization options
options = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations',1e5);

% Call fminimax
[PV_opt, fval] = lsqnonlin(objective, PV_axis_0, [], [], [], [], [], [], nonlcon, options);

% Extract optimized P_axis and V_axis
P_axis_opt = reshape(PV_opt(1:3*num_axis), [3, num_axis]);
V_axis_opt = reshape(PV_opt(3*num_axis+1:end), [3, num_axis]);

residual_after = residual_fun(P_axis_opt, V_axis_opt, magnet_conf(1:3,:), magnet_conf(4:6,:), cylinder_volume, B_r, B_star);

%%
% mag_field(P_axis_0, V_axis_0, [0;0;0.1],[0;0;1],B_r,cylinder_volume)*1e6
mag_field(P_axis_0, V_axis_0, [0;0;0.1],[0;-0.7071;0.7071],B_r,cylinder_volume)*1e6

%% Residual function
function residual = residual_fun(P_axis, V_axis, P_magnet, Mu_hat_magnet, Volumn, B_r, B_star)
    B_star = reshape(B_star, [], 1);
    B_collection = mag_field_vector(P_axis, V_axis, P_magnet, Mu_hat_magnet, Volumn, B_r);
    B_collection = reshape(B_collection, [], 1) * 1e6;

    residual = abs(B_star - B_collection);
end

%% Constraint function
function [c, ceq] = constraints(PV, ncols)
    % Extract P_axis and V_axis from PV
    P_axis = reshape(PV(1:3*ncols), [3, ncols]);
    V_axis = reshape(PV(3*ncols+1:end), [3, ncols]);

    % Equality constraints for V_axis norms
    ceq = zeros(1, ncols);
    for i = 1:ncols
        ceq(i) = norm(V_axis(:, i)) - 1;
    end

    % No inequality constraints
    c = [];
end

%% Functions
function B_collection = mag_field_vector(P_axis, V_axis, P_magnet, Mu_hat_magnet, Volumn, B_r)
    % P_magnet: positions of the magnet for n instances (3*n)
    % Mu_magnet: orientations of the magnet for n instances (3*n)
    % P_axis: positions of the axis for axis number m (3*m)
    % V_axis: orientations of the axis for axis number m (3*m)
    % Volumn: volumn of the magnet
    % B_r: magnetic flux density of the magnet

    % B_collection: for the timeframe of n instances for
    % P_magnet and Mu_hat_magnet, magnetic field for each sensor axis (3m*n)
    
    B_collection = [];
    % For magnet configuration, construct the magnetic field reading
    for i = 1:size(P_magnet,2)
        % Extract the magnet position and orientation
        p_magnet = P_magnet(:,i);
        mu_hat_magnet = Mu_hat_magnet(:,i);
        
        % Construct the sensor reading for this time stamp's magnet
        % position and orientation
        B_collection = [B_collection, mag_field(P_axis, V_axis, p_magnet, mu_hat_magnet, B_r, Volumn)];
    end
end

function B_vector = mag_field(P_axis, V_axis, p_magnet, mu_hat_magnet, B_r, Volumn)
    % p_magnet: magnet's position (3*1)
    % mu_hat_magnet: magnet's orientation (3*1, unit vector)
    % P_axis: sensor axis position (3*3m)
    % V_axis: sensor axis orientation (3*3m)

    B_vector = [];
    % For each axis position and orientation, construct the magnetic field
    for i = 1:size(P_axis, 2)
        p_axis = P_axis(:,i);
        v_axis = V_axis(:,i);
        
        % Magnetic field for one axis
        B_axis = mag_field_axis(p_axis, v_axis, p_magnet, mu_hat_magnet, B_r, Volumn);
        
        % Put that into the collection
        B_vector = [B_vector; B_axis];        
    end
end

function B_axis = mag_field_axis(p_axis, v_axis, p_magnet, mu_hat_magnet, B_r, Volumn)
    r = p_axis-p_magnet;
    r_hat = r/norm(r);
    B = (B_r*Volumn)/(4*pi*(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*mu_hat_magnet;

    % Project the field onto the sensor's axis
    B_axis = dot(B, v_axis);
end