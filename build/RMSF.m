function [RMS] = RMSF(x,sur,p)
    
    if (nargin==2)
        p = sur;
        sur=false;
    end
    
    %   Scaling the function inputs from unity 
    
    lb = p{2};
    ub = p{3};
    n_postrack = p{4};
    n_memory = p{5};
    
    y=scaling(x',lb,ub,2);
        
    c = y(1);
    lambda = y(2);
    
    %   Calling the Attitude estimator
    if ~(sur)
        RMS = Main_delay_fun(lambda, c, n_postrack, n_memory);
    else
        [RMS,std,ei,cdf] = sgtelib_server_predict(y);
    end
 
    
end