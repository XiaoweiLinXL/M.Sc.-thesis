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
folder = fullfile(pwd, "data_10cm");

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
        B_star = [B_star, mean(data.user_data.Field_filt(:, data.save_element(j)-10: data.save_element(j)), 2)];

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