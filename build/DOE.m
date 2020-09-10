clc
clear all
close all
format compact
format long
addpath('Support_functions')

p = 2000;  % Number of points
N = 2;   % Number of dimensions

% DOE 1
ll_c = 8; ul_c = 20;
ll_lambda = 0.01; ul_lambda = 5.0;
% DOE 2
% ll_c = 0.0; ul_c = 20;
% ll_lambda = 0.0; ul_lambda = 2.5;
% DOE 3
% ll_c = 5.0; ul_c = 15.0;
% ll_lambda = 0.01; ul_lambda = 5.0;

lb = [ll_c ll_lambda]; % lower bounds for lambda , C
ub = [ul_c ul_lambda]; % upper bounds for lambda , C

n_postrack = 20; n_memory = 30;
param = {1,lb,ub,n_postrack,n_memory};

X = lhsdesign(p,N,'criterion','correlation'); % generate latin hypercube samples

for i=1:p 
    [Err(i),v_Err(i),Fcont(i)] = RMSF_cpp(X(i,:)',param);
    i
end

%% Save the DOE results
save('DOE_results_V1','X','Err','v_Err','lb','ub')