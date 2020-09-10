clearvars
close all
clc
addpath SGTE_matlab_server
addpath Support_functions
addpath ./Support_functions/hatchfill2_r8
addpath ./Support_functions/export_fig

%% Problem definition
load('DOE_results_V3.mat','X','Err','v_Err','lb','ub');

lhs_data_normalize = X;
obj_data = v_Err;

%% Design space sampling
n_samples = 25; % <-------------------------------------------------------- SET NUMBER OF SAMPLES FOR DOE
lhs_data = scaling(lhs_data_normalize, lb, ub, 2);
%% Construct surrogate models <---------------------------------------------------- CHOOSE DIFFERENT SURROGATE MODELS
%-------------------------------------------------------------------------%
% For default hyperparameters
% model = 'TYPE LOWESS DEGREE 2 KERNEL_TYPE D1 KERNEL_SHAPE 1.12073 DISTANCE_TYPE NORM2 RIDGE 0.0125395';
% model = 'TYPE KRIGING RIDGE 1.01723e-16 DISTANCE_TYPE NORM2 METRIC OECV BUDGET 200';
% For Hyperparameter Optimization
budget = 200; out_file = 'surrogate_model.sgt';
% model = ['TYPE LOWESS ', 'DEGREE OPTIM ', 'RIDGE OPTIM ', 'KERNEL_TYPE OPTIM ', 'KERNEL_COEF OPTIM ', 'DISTANCE_TYPE OPTIM ', 'METRIC OECV ', 'BUDGET ', num2str(budget), ' OUTPUT ', out_file];
% model = ['TYPE KS ', 'KERNEL_TYPE OPTIM ', 'KERNEL_COEF OPTIM ', 'DISTANCE_TYPE OPTIM ', 'METRIC OECV ','BUDGET ', num2str(budget), ' OUTPUT ', out_file];
% model = ['TYPE RBF ', 'KERNEL_TYPE OPTIM ', 'KERNEL_COEF OPTIM ', 'DISTANCE_TYPE OPTIM ', 'RIDGE OPTIM ', 'METRIC OECV ', 'BUDGET ', num2str(budget), ' OUTPUT ', out_file];
% model = ['TYPE KRIGING ', 'RIDGE OPTIM ', 'DISTANCE_TYPE OPTIM ', 'METRIC OECV ', 'BUDGET ', num2str(budget), ' OUTPUT ', out_file];
model = ['TYPE ENSEMBLE ', 'WEIGHT OPTIM ', 'METRIC OECV ', 'DISTANCE_TYPE OPTIM ','BUDGET ', num2str(budget),' OUTPUT ', out_file];
%-------------------------------------------------------------------------%

sgtelib_server_start(model,true,true)
% Test if server is ok and ready
sgtelib_server_ping;
% Feed server
sgtelib_server_newdata(lhs_data,obj_data');

% metric_str = 'OECV';
% metric = sgtelib_server_metric(metric_str);
% fprintf('===============================\n')
% fprintf('The OECV is : %f\n',metric)
% metric_str = 'RMSECV';
% metric = sgtelib_server_metric(metric_str);
% fprintf('The RMSECV is : %f\n',metric)
% fprintf('===============================\n')

%Prediction

% lb_plot = [13 1.6];
% ub_plot = [21.5 2.3];
% lb_plot = [9.0 0.1];
% ub_plot = [20.0 5.0];
lb_plot = lb;
ub_plot = ub;

res = 70;
X = gridsamp([lb_plot; ub_plot], res);
[YX,std,ei,cdf] = sgtelib_server_predict(X);
% sgtelib_server_stop; %stop the server 

X1 = reshape(X(:,1),res,res); X2 = reshape(X(:,2),res,res);
YX = reshape(YX, size(X1));

%% Plot design space

fig1 = figure(1);
h = axes(fig1);
axis(h,[lb_plot(1),ub_plot(1),lb_plot(2),ub_plot(2)]) % fix the axis limits

[cc, hh] = contourf(h,X1, X2, YX,20); % plot contour
% hold on
% plot(lhs_data(:,1),lhs_data(:,2),'.k','markersize',10)
colorbar(h)
xlabel('relay amplitude (c)','interpreter','latex','fontsize',16)
ylabel('($\lambda$)','interpreter','latex','fontsize',16)
set(fig1,'color','w');
export_fig('surrogate_function.pdf','-p0.002',fig1); 
export_fig('surrogate_function.png','-p0.002','-r600',fig1); 