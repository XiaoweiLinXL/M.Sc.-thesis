close all
clear 
clc

%% Constants
sphere_diameter = 6.35e-3; % 6.35mm, diameter of the spherical magnet
sphere_volumn = (4/3)*pi*(sphere_diameter/2)^3; % Volumn of a sphere 

B_r = 1.32; % Magnetic flux density (T)
num_sphere = 10; % Number of spherical magnets
Volumn = sphere_volumn * num_sphere; % Total volumn of the magnets

%%
load("data.mat")

%%
p_magnet = P0(:,1:13);
p_magnet(1,:) = p_magnet(1,:) - 5*sphere_diameter;
mu_hat_magnet = repmat([-1;0;0], 1, size(p_magnet, 2));

field = [];
for i = 1:size(p_magnet, 2)
    field = [field, mag_field_timestamp(p_magnet(:,i), mu_hat_magnet(:,i), sens_pos, sens_or, B_r, Volumn)];
end

field = field * 1e6;

measured_field = user_data.Field(:,save_element(1:13),:);

figure
tiledlayout(2,2)
for i = 1:4
    nexttile
    error = field(3*i-2:3*i,:) - measured_field(:,:,i);
    plot(1:size(field,2), error)
end

%% Functions
function B = mag_field_vector_timeframe(P_magnet, Mu_hat_magnet, P_axis, V_axis, Volumn, B_r)
    % P_magnet: positions of the magnet for time stamps n (3*n)
    % Mu_magnet: orientations of the magnet for time stamps n (3*n)
    % P_axis: positions of the axis for axis number m (3*m)
    % V_axis: orientations of the axis for axis number m (3*m)
    % Volumn: volumn of the magnet
    % B_r: magnetic flux density of the magnet
    
    % For each time stamp, construct the magnetic field reading
    for i = 1:size(P_magnet,2)
        % Extract the magnet position and orientation
        p_magnet = P_magnet(:,i);
        mu_hat_magnet = Mu_hat_magnet(:,i);
        
        % For each sensor axis, construct the field reading


        

    end

end

function B_timestamp = mag_field_timestamp(p_magnet, mu_hat_magnet, P_axis, V_axis, B_r, Volumn)
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