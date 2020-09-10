function [RMS,v_RMS,FCONT] = RMSF_cpp(x,sur,p)
    
    if (nargin==2)
        p = sur;
        sur=false;
    end
    
    %   Scaling the function inputs from unity 
    
    lb = p{2};
    ub = p{3};
    n_postrack = p{4};
    n_memory = p{5};
    
    X_unscaled = scaling(x',lb,ub,2);
    c = X_unscaled(1);
    lambda = X_unscaled(2);
    
    %   Calling the Attitude estimator
    if ~(sur)
        command = ['cd',' ',pwd,'&',' ','IMU_Attitude_Estimator_Ex.exe', ' ', num2str(c), ' ', num2str(lambda), ' ', num2str(n_postrack), ' ', num2str(n_memory)];
        status = system(command); 

        if status == 0
            
            % RMS error in position metric
            fileID = fopen('POS_ERR_RMS.txt','r');
            RMS = fscanf(fileID,'%f');
            fclose(fileID);
            
            % RMS error in velocity metric
            fileID = fopen('VEL_ERR_RMS.txt','r');
            v_RMS = fscanf(fileID,'%f');
            fclose(fileID);
            
            % custom error metric
            fileID = fopen('F_CONT_ABS.txt','r');
            FCONT = fscanf(fileID,'%f');
            fclose(fileID);
            else
            error('Application did not execute succesfully');
        end
    else
        [Y,std,ei,cdf] = sgtelib_server_predict(y);
        [RMS, v_RMS, FCONT] = Y;
    end
 
    
end