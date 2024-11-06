close all
clear 
clc

%% Constants
sphere_diameter = 6.35e-3; % 6.35mm, diameter of the spherical magnet
sphere_volumn = (4/3)*pi*(sphere_diameter/2)^3; % Volumn of a sphere 

B_r = 1.32; % Magnetic flux density (T)
num_sphere = 10; % Number of spherical magnets
Volumn = sphere_volumn * num_sphere; % Total volumn of the magnets

num_timestamp = 13;

%%
load("data.mat")

%%
P_magnet_total = P0; % ball chain positions at time stamps save_element
P_magnet_total(1,:) = P_magnet_total(1,:) - 5*sphere_diameter; % offset the dipole to the middle of the chain 
Mu_hat_magnet_total = repmat([-1;0;0], 1, size(P_magnet_total, 2)); % dipole directions at save_element are facing -x
B_star_total = user_data.Field(:,save_element,:);
B_star_total = [B_star_total(:,:,1);B_star_total(:,:,2);B_star_total(:,:,3);B_star_total(:,:,4)];

% Select the time frame we want to use for calibration
num_timeframe = size(P_magnet_total, 2) - num_timestamp;

for i = 1:num_timeframe
    % Magnet position for time stamps in the selected time frame
    P_magnet = P_magnet_total(:,num_timestamp*i - (num_timestamp-1):num_timestamp*i);
    % Magnet orientation for time stamps in the selected time frame
    Mu_hat_magnet = Mu_hat_magnet_total(:,num_timestamp*i - (num_timestamp-1):num_timestamp*i);
    % Measured magnetic field
    B_star = B_star_total(:,num_timestamp*i - (num_timestamp-1):num_timestamp*i);


    residual_fun(sens_pos, sens_or, P_magnet, Mu_hat_magnet, Volumn, B_r, B_star)
end



%%
function residual = residual_fun(P_axis, V_axis, P_magnet, Mu_hat_magnet, Volumn, B_r, B_star)
    B_star = reshape(B_star, [], 1);
    B_timeframe_collection = mag_field_vector_timeframe(P_axis, V_axis, P_magnet, Mu_hat_magnet, Volumn, B_r);
    B_timeframe_collection = reshape(B_timeframe_collection, [], 1) * 1e6;

    residual = B_star - B_timeframe_collection;
end


%% Functions
function B_timeframe_collection = mag_field_vector_timeframe(P_axis, V_axis, P_magnet, Mu_hat_magnet, Volumn, B_r)
    % P_magnet: positions of the magnet for time stamps n (3*n)
    % Mu_magnet: orientations of the magnet for time stamps n (3*n)
    % P_axis: positions of the axis for axis number m (3*m)
    % V_axis: orientations of the axis for axis number m (3*m)
    % Volumn: volumn of the magnet
    % B_r: magnetic flux density of the magnet

    % B_timeframe_collection: for the timeframe of n timestamps for
    % P_magnet and Mu_hat_magnet, magnetic field for each sensor axis (3m*n)
    
    B_timeframe_collection = [];
    % For each time stamp, construct the magnetic field reading
    for i = 1:size(P_magnet,2)
        % Extract the magnet position and orientation
        p_magnet = P_magnet(:,i);
        mu_hat_magnet = Mu_hat_magnet(:,i);
        
        % Construct the sensor reading for this time stamp's magnet
        % position and orientation
        B_timeframe_collection = [B_timeframe_collection, mag_field_timestamp(P_axis, V_axis, p_magnet, mu_hat_magnet, B_r, Volumn)];
    end
end

function B_timestamp = mag_field_timestamp(P_axis, V_axis, p_magnet, mu_hat_magnet, B_r, Volumn)
    % p_magnet: magnet's position (3*1)
    % mu_hat_magnet: magnet's orientation (3*1, unit vector)
    % P_axis: sensor axis position (3*3m)
    % V_axis: sensor axis orientation (3*3m)

    B_timestamp = [];
    % For each axis position and orientation, construct the magnetic field
    for i = 1:size(P_axis, 2)
        p_axis = P_axis(:,i);
        v_axis = V_axis(:,i);
        
        % Magnetic field for one axis
        B_axis = mag_field_axis(p_axis, v_axis, p_magnet, mu_hat_magnet, B_r, Volumn);
        
        % Put that into the collection
        B_timestamp = [B_timestamp; B_axis];        
    end
end

function B_axis = mag_field_axis(p_axis, v_axis, p_magnet, mu_hat_magnet, B_r, Volumn)
    r = p_axis-p_magnet;
    r_hat = r/norm(r);
    B = (B_r*Volumn)/(4*pi*(norm(r)^3))*(3*(r_hat*r_hat.')-eye(3))*mu_hat_magnet;

    % Project the field onto the sensor's axis
    B_axis = dot(B, v_axis);
end