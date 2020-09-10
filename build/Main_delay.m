
%% Start of script
addpath('Support_functions')
addpath('quaternion_library');      % include quaternion library
addpath('datasets');      % include datasets
close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal

%% Read Data File
count=1;

fid = 'NAV3_data.bin';
data = load(fid, '-ascii');
% data = data(1:1:1000,:);
Rows = size(data,1);
measurements = data(:,:);
time_raw =  data(:,25);
groundtruth_ds = data(:,[10 11 12]);
groundtruth_pos_ds = data(:,[13 14 15]);
groundtruth_true = data(:,[16 17 18]);
%%groundtruth_true = groundtruth_true()* 180 / pi;
groundtruth_pos_true = data(:,[19 20 21]);
groundtruth_vel_true = data(:,[22 23 24]);
acceleration = data(:,[1 2 3]);
gyroscope = data(:,[4 5 6]);
magnetometer = data(:,[7 8 9]);
quaternion = zeros(length(time_raw), 4);

%% Variable Initialization
Kp = 2.0;
Ki = 1.0;
c = 15; % as per the optimization c=15.09 lambda1.914
lambda = 2;
freq = 1000.0;
g = 9.81;
p = 2.0; 
q= 1.0;
dt = 0.001;
pos_err_sq_sum = 0.0;
pos_err_RMS = 0;
integ_vel = zeros(1,3);
integ_pos = zeros(1,3);
% integ_vel = zeros();
% integ_pos = zeros();

gravity_i = [0.0; 0.0; g];
n_memory = 21;
n_postrack = 40;
%% Process sensor data through algorithm


%AHRS = MadgwickAHRS('SamplePeriod',1/freq, 'Kp', Kp,'Ki',Ki);
AHRS = MahonyAHRS('SamplePeriod',1/freq, 'Kp', Kp,'Ki',Ki);
i=1;

% initialize pos vector 
pos_est_current = groundtruth_pos_ds(i, :); 

% initialize acc vector
acc_data = acceleration(1,:);
gyro_data = gyroscope(1,:);
[acc_i(1,:),euler(1,:)] = compute_acc_i(AHRS,gyro_data,acc_data,gravity_i);

acc_est(i, :) = acc_i(i, :);
acc_est_prev = acc_est(i, :);

% initialize super_twist_vel vector
vel_inc = acc_est(i, :);
super_twist_vel(i, :) = (vel_inc * dt);
st_vel_prev = super_twist_vel(i, :);  % get last row of the matrix


while(i<Rows+1)
    
    %% Mahony section
    acc_data = acceleration(i,:);
    gyro_data = gyroscope(i,:);
    [acc_i(i,:),euler(i,:)] = compute_acc_i(AHRS,gyro_data,acc_data,gravity_i);
    
    %% CH Boiko section
    
    acc_i_raw = acc_i(i, :);
    pos_true = groundtruth_pos_ds(i, :);

    [acc_est(i, :), vel_est(i, :), super_twist_vel(i, :), pos_est(i, :), integ_vel, integ_pos] = chib_filter(...
        pos_est_current, pos_true, acc_i_raw, acc_est_prev, st_vel_prev, integ_vel, integ_pos, c, dt, lambda, q, p);
    
    if rem(i,n_postrack) == 1 && i > n_postrack
        
        
        if i - n_memory == 0 % First correction block
            j = 1;
        else
            j = i-n_memory;
        end
        
        % Reset the integrators
        integ_vel = vel_est(j,:);
        integ_pos = pos_est(j,:);
        %% initialize for next loop iteration 
        % pos vector
        pos_est_current = pos_est(j, :); % get last row of the matrix
        % pos_est will be updated at end of loop!
        % acc vector
        acc_est_prev = acc_est(j,:); 
        % super_twist_vel vector
        st_vel_prev = super_twist_vel(j, :);  % get last row of the matrix
        
        while (j <= i)

            acc_i_raw = acc_i(j, :);

            [acc_est(j, :), vel_est(j, :), super_twist_vel(j, :), pos_est(j, :), integ_vel, integ_pos] = chib_filter(...
            pos_est_current, pos_true, acc_i_raw, acc_est_prev, st_vel_prev, integ_vel, integ_pos, c, dt, lambda, q, p);
            
            % fprintf('x, y, z: %i | %i | %.14f | %.14f | %.14f\n',i,j,pos_true(1),pos_true(2),pos_true(3))
        
            j = j + 1;
        
            %% initialize for next loop iteration 
            % pos vector
            pos_est_current = pos_est(j - 1, :); % get last row of the matrix
            % pos_est will be updated at end of loop!
            % acc vector
            acc_est_prev = acc_est(j - 1,:); 
            % super_twist_vel vector
            st_vel_prev = super_twist_vel(j - 1, :);  % get last row of the matrix
            
        end
    end
    
    i=i+1;
    
    %% initialize for next loop iteration 
    % pos vector
    pos_est_current = pos_est(i - 1, :); % get last row of the matrix
    % pos_est will be updated at end of loop!
    % acc vector
    acc_est_prev = acc_est(i - 1,:); 
    % super_twist_vel vector
    st_vel_prev = super_twist_vel(i - 1, :);  % get last row of the matrix
       
end

%% Compute RMS error with respect to true position

pos_err_true = pos_est - groundtruth_pos_true;
pos_err_sq = pos_err_true.^2;
pos_err_sq_sum = sum(pos_err_sq,2);
pos_err_sq_sum_sqrt = sqrt(pos_err_sq_sum);
pos_err_RMS = sum(pos_err_sq_sum_sqrt)/i;

% pos_err_sq_sum_sqrt = 0;
% i = 1;
% while(i<Rows+1)
% %% Compute RMS error with respect to true position
% 
% pos_err_true = pos_est(i,:) - groundtruth_pos_true(i,:);
% pos_err_sq = pos_err_true.^2;
% pos_err_sq_sum = sum(pos_err_sq);
% pos_err_sq_sum_sqrt = pos_err_sq_sum_sqrt + sqrt(pos_err_sq_sum);
% 
% i = i + 1;
% 
% end
% pos_err_RMS = sum(pos_err_sq_sum_sqrt)/i;

%% Plot algorithm output position

figure('Name', 'Position');
axis(1) = subplot(1,3,1);
hold on;
plot(  time_raw, pos_est(:,1), 'r');
plot(  time_raw, groundtruth_pos_ds(:,1), 'g');
plot(  time_raw, groundtruth_pos_true(:,1), 'k');
legend('Est.', 'DownSampled','True');
xlabel('Time (s)');
ylabel('Meter');
title('X');
hold off;

axis(2) = subplot(1,3,2);
hold on;
plot(  time_raw, pos_est(:,2), 'r');
plot(  time_raw, groundtruth_pos_ds(:,2), 'g');
plot(  time_raw, groundtruth_pos_true(:,2),'k');
legend('Est.', 'DownSampled','True');
xlabel('Time (s)');
ylabel('Meter');
title('Y');
hold off;

axis(3) = subplot(1,3,3);
hold on;
plot(  time_raw, pos_est(:,3), 'r');
plot(  time_raw, groundtruth_pos_ds(:,3), 'g');
plot(  time_raw, groundtruth_pos_true(:,3), 'k');
legend('Est.', 'DownSampled','True');
xlabel('Time (s)');
ylabel('Meter');
title('Z');
hold off;
linkaxes(axis, 'x');

% export estimation results
save('pos_est.mat','pos_est','time_raw','groundtruth_pos_ds','groundtruth_pos_true','acc_i','euler')

%% Plot algorithm output velocity
% 
% figure('Name', 'Velocity');
% axis(1) = subplot(1,3,1);
% hold on;
% plot(time_raw, super_twist_vel(:, 1), 'r');
% legend('Est.');
% xlabel('Time (s)');
% ylabel('m/s');
% title('X');
% hold off;
% 
% axis(2) = subplot(1,3,2);
% hold on;
% plot(time_raw, super_twist_vel(:, 2), 'r');
% legend('Est.');
% xlabel('Time (s)');
% ylabel('m/s');
% title('Y');
% hold off;
% 
% axis(3) = subplot(1,3,3);
% hold on;
% plot(time_raw, super_twist_vel(:, 3), 'r');
% legend('Est.');
% xlabel('Time (s)');
% ylabel('m/s');
% title('Z');
% hold off;
% linkaxes(axis, 'x');
%% Plot algorithm output as Euler angles
% The first and third Euler angles in the sequence (phi and psi) become
% unreliable when the middle angles of the sequence (theta) approaches ±90
% degrees. This problem commonly referred to as Gimbal Lock.
% See: http://en.wikipedia.org/wiki/Gimbal_lock

figure('Name', 'Euler Angles');
axis(1) = subplot(1,3,1);
hold on;
plot(time_raw, euler(:,1), 'r');
plot(time_raw, (groundtruth_true(:,1)*180/pi), 'g');
legend('Est.', 'True');
xlabel('Time (s)');
ylabel('Angle (deg)');
title('\phi');
hold off;

axis(2) = subplot(1,3,2);
hold on;
plot(time_raw, euler(:,2), 'k');
plot(time_raw, (groundtruth_true(:,2)*180/pi), 'g');
legend('Est.', 'True');
xlabel('Time (s)');
ylabel('Angle (deg)');
title('\theta');
hold off;

axis(3) = subplot(1,3,3);
hold on;
plot(time_raw, euler(:,3), 'b');
plot(time_raw, (groundtruth_true(:,3)*180/pi), 'g');
legend('Est.', 'True');
xlabel('Time (s)');
ylabel('Angle (deg)');
title('\psi');
hold off;
linkaxes(axis, 'x');
%% End of script

function [acc_i_out,euler_out] = compute_acc_i(AHRS,gyro_data,acc_data,gravity_i)
    
    g = 9.81;

    AHRS.UpdateIMU(gyro_data, acc_data);	% gyroscope units must be radians & UpdateIMU for system without paraometer / Update for system with paraometers
    quaternion = AHRS.Quaternion;
    euler = quatern2euler(quaternConj(quaternion)) * (180/pi);	% use conjugate for sensor frame relative to Earth and convert to degrees.
    Rotation_matrix = quatern2rotMat(quaternion);
        
    gravity_b = ( Rotation_matrix ) * gravity_i;

    acc_b = (acc_data*g) - gravity_b';
    acc_i = Rotation_matrix \ acc_b';
    
    euler_out = euler;
    acc_i_out = acc_i;
    
end

function [acc_est_out, vel_est_out, super_twist_vel_out, pos_est_out, integ_vel, integ_pos] = chib_filter(...
    pos_est_current, pos_true, acc_i_raw, acc_est_prev, st_vel_prev, integ_vel, integ_pos, c, dt, lambda, q, p)

    %% Acceleration calculations
    pos_err = pos_est_current - pos_true;
    pos_err_neg = pos_err * -1;

    relay = sign(pos_err_neg); 
    
    acc_est = acc_i_raw + (c * relay);
    
    acc_est_current = acc_est;
    %% Velocity calculations

    vel_inc = (acc_est_prev + acc_est_current) / 2;
    integ_vel = (vel_inc * dt) + integ_vel; %% integrate acceleration to get velocity
    vel_est = integ_vel;

    %% Lipshitchz discontinous function
    super_twist = abs(pos_err_neg) .^ (q / p) .* sign( pos_err_neg);
    %Linear function
%     super_twist = pos_err_neg;

    super_twist_vel = vel_est + (lambda * super_twist);
    st_vel_current = super_twist_vel;
    %% Position calculations

    pos_inc = (st_vel_prev + st_vel_current) / 2;
    integ_pos = (pos_inc * dt) + integ_pos; %% integrate acceleration to get velocity
    pos_est = integ_pos;
    
	acc_est_out = acc_est;
    vel_est_out = vel_est;
	super_twist_vel_out = super_twist_vel;
    pos_est_out = pos_est;

end
