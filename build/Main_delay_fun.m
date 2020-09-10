function pos_err_RMS = Main_delay_fun(lambda, c, n_postrack, n_memory)

    %% Start of script

    addpath('quaternion_library');      % include quaternion library
    addpath('datasets');      % include datasets

    %% Read Data File
    count=1;

    fid = 'NAV3_data.bin';
    data = load(fid, '-ascii');
    % data = data(1:1:1000,:);
    Rows = size(data,1);
    measurements = data(:,:);
    time_raw =  data(:,22);
    groundtruth_ds = data(:,[10 11 12]);
    groundtruth_pos_ds = data(:,[13 14 15]);
    groundtruth_true = data(:,[16 17 18]);
    %%groundtruth_true = groundtruth_true()* 180 / pi;
    groundtruth_pos_true = data(:,[19 20 21]);
    acceleration = data(:,[1 2 3]);
    gyroscope = data(:,[4 5 6]);
    magnetometer = data(:,[7 8 9]);
    quaternion = zeros(length(time_raw), 4);

    %% Variable Initialization
    Kp = 2.0;
    Ki = 1.0;
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
                pos_true = groundtruth_pos_ds(i, :);

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

    %% End of script

end

function [acc_i_out,euler_out] = compute_acc_i(AHRS,gyro_data,acc_data,gravity_i)
    
    g = 9.81;

    AHRS.UpdateIMU(gyro_data, acc_data);	% gyroscope units must be radians & UpdateIMU for system without paraometer / Update for system with paraometers
    quaternion = AHRS.Quaternion;
    euler = quatern2euler(quaternConj(quaternion)) * (180/pi);	% use conjugate for sensor frame relative to Earth and convert to degrees.
    Rotation_matrix = quatern2rotMat(quaternion);
        
    gravity_b = Rotation_matrix * gravity_i;

    acc_b = acc_data * g - gravity_b';
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
    super_twist =abs(pos_err_neg) .^ (q / p) .* sign( pos_err_neg);
    %Linear function
    %super_twist = pos_err_neg;
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
