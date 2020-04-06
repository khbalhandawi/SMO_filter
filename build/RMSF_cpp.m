function [RMS,FCONT] = RMSF_cpp(x,sur,p)
    
    if (nargin==2)
        p = sur;
        sur=false;
    end
    
    %   Scaling the function inputs from unity 
    
    n_postrack = p{2};
    n_memory = p{3};
    
    c = x(1);
    lambda = x(2);
    
    %   Calling the Attitude estimator
    if ~(sur)
        command = ['cd',' ',pwd,'&',' ','IMU_Attitude_Estimator_Ex.exe', ' ', num2str(c), ' ', num2str(lambda), ' ', num2str(n_postrack), ' ', num2str(n_memory)];
        status = system(command); 

        if status == 0
            % fileID = fopen('F_CONT_ABS.txt','r');
            fileID = fopen('POS_ERR_RMS.txt','r');
            RMS = fscanf(fileID,'%f');
            fclose(fileID);
            fileID = fopen('F_CONT_ABS.txt','r');
            FCONT = fscanf(fileID,'%f');
            fclose(fileID);
            else
            error('Application did not execute succesfully');
        end
    else
        [RMS,std,ei,cdf] = sgtelib_server_predict(y);
    end
 
    
end

function [RMS] = RMSF(x,p)
    
    %   Scaling the function inputs from unity 
    lb = p{1};
    ub = p{2};
    y=scaling(x,lb,ub,2);
    c = y(1);
    lambda = y(2);
    

    
end