%% Start of script
addpath('Support_functions')
addpath('Support_functions\export_fig')
addpath('quaternion_library');      % include quaternion library
addpath('datasets');      % include datasets
close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal

%% Call cpp application
% lambda = 2.039450;
% c = 5.324800;

lambda = 10.0;
c = 7.0;

% c = 1.44035;
% lambda = 1.25285;

% lambda = 0.2174;
% c = 6.377;

ll_c = 0.0; ul_c = 20;
ll_lambda = 0.0; ul_lambda = 5.0;
lb = [ll_c ll_lambda]; % lower bounds for lambda , C
ub = [ul_c ul_lambda]; % upper bounds for lambda , C

n_delay = 20;
n_memory = 30 ; % how many anti-delay units do you want
n_postrack = 0;

X = [c, lambda];
X = scaling(X',lb,ub,1);

param = {1,lb,ub,n_postrack,n_memory};
[Err,v_Err,Font] = RMSF_cpp(X',param);

export_fig = false;

%% Read Data File
count=1;
% range = [31750:32050]; % static
% range = [28600:29842]; % dynamic

fid = 'time_raw.bin';
data = load(fid, '-ascii');
range = [1:length(data(:,1))]; % dynamic
time_raw = data(range,1);

fid = 'pos_est.bin';
data = load(fid, '-ascii');
Rows = size(data,1);
pos_est_cpp = data(range,[1 2 3]);

fid = 'acc_est.bin';
data = load(fid, '-ascii');
acc_i_cpp = data(range,[1 2 3]);

fid = 'vel_est.bin';
data = load(fid, '-ascii');
vel_est_cpp = data(range,[1 2 3]);

fid = 'NAV3_data.bin';
data = load(fid, '-ascii');
groundtruth_pos_ds = data(range,[13 14 15]);
groundtruth_pos_true = data(range,[19 20 21]);
groundtruth_vel_true = data(range,[22 23 24]);
groundtruth_acc_true = data(range,[1 2 3]);

% %% Plot algorithm output position
% 
% fig1 = figure('Name', 'Position');
% axis(1) = subplot(1,3,1);
% hold on;
% plot(time_raw, vel_est_cpp(:,1), '-r', 'linewidth',1.5);
% % plot(time_raw, groundtruth_pos_ds(:,1), '-g', 'linewidth',1.5);
% plot(time_raw, groundtruth_vel_true(:,1), '-k', 'linewidth',1.5);
% xlabel('Time (s)','interpreter','latex');
% ylabel('m/s','interpreter','latex');
% title('$v_x$ (m/s)','interpreter','latex');
% % ylim([1.08,1.13])
% hold off;
% 
% axis(2) = subplot(1,3,2);
% hold on;
% plot(time_raw, vel_est_cpp(:,2), '-r', 'linewidth',1.5);
% % plot(time_raw, groundtruth_pos_ds(:,2), '-g', 'linewidth',1.5);
% plot(time_raw, groundtruth_vel_true(:,2),'-k', 'linewidth',1.5);
% xlabel('Time (s)','interpreter','latex');
% ylabel('m/s','interpreter','latex');
% title('$v_y$ (m/s)','interpreter','latex');
% % ylim([-1.548,-1.536])
% % ylim([-1.17196,-0.3343])
% hold off;
% 
% axis(3) = subplot(1,3,3);
% hold on;
% est_h = plot(time_raw, vel_est_cpp(:,3), '-r', 'linewidth',1.5);
% % ds_h = plot(time_raw, groundtruth_pos_ds(:,3), '-g', 'linewidth',1.5);
% true_h = plot(time_raw, groundtruth_vel_true(:,3), '-k', 'linewidth',1.5);
% xlabel('Time (s)','interpreter','latex');
% ylabel('m/s','interpreter','latex');
% title('$v_z$ (m/s)','interpreter','latex');
% % ylim([1.13,1.19])
% hold off;
% linkaxes(axis, 'x');
% 
% %=========================== FIGURE SETTINGS =============================%
% %% Figure resizing
% ax_h = -0.08; ax_bot = 0; ann_pos = 0.45;    %<-------- Edit as necessary to fit figure properly
% fig_width = 1200; fig_height = 500;
% ch = fig1.Children;
% 
% set(fig1, 'Position', [50, 50, fig_width, fig_height])
% for par = 1:1:3
%     
%     % Plot points
%     ax = ch(par);
%     
%     x = get(ax,'XTickLabel');
%     set(ax,'XTickLabel',x,'FontName','Times','fontsize',14)
%     y = get(ax,'XTickLabel');
%     set(ax,'XTickLabel',y,'FontName','Times','fontsize',14)
%     
%     outerpos = ax.OuterPosition;
% 
%     ti = ax.TightInset;
%     ti(1:2:3) = ti(1:2:3) - 0.015;
%     ti(2:2:4) = ti(2:2:4) + 0.058;
% 
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4) + ax_h;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2) - ax_bot;
% 
%     pos = [left bottom ax_width ax_height];
%     ax.Position = pos;
%     
% end
% 
% legend_items = [est_h true_h];
% legend_entries = [{'Estimated Velocity'}, {'True Velocity'}];
% 
% lh = legend(legend_items,legend_entries,'Orientation','horizontal');
% % lh = legend([h_obj lc_h],[title_obj label_lc],'Orientation','horizontal');
% rect = [0.3813, 0.92, .27, .0528];
% set(lh, 'Position', rect, 'interpreter', 'latex', 'fontsize', 16)
% 
% if export_fig
%     set(fig1,'color','w');
%     savefig(fig1,'./Sample_plots/fig1.fig')
%     export_fig(['./Sample_plots/response_vel_delay_',num2str(n_delay),'_c_',num2str(c),'_l_',num2str(lambda),'.pdf'],'-p0.002',fig1); 
%     export_fig(['./Sample_plots/response_vel_delay_',num2str(n_delay),'_c_',num2str(c),'_l_',num2str(lambda),'.png'],'-p0.002','-r600',fig1); 
% end
%% Plot algorithm output position

fig2 = figure('Name', 'Position');
hold on;
est_h = plot(time_raw, pos_est_cpp(:,2), '-r', 'linewidth',1.5);
ds_h = plot(time_raw, groundtruth_pos_ds(:,2), '-g', 'linewidth',1.5);
true_h = plot(time_raw, groundtruth_pos_true(:,2), '-k', 'linewidth',1.5);
xlabel('Time (s)','interpreter','latex');
ylabel('y (m)','interpreter','latex');
xlim([min(time_raw),max(time_raw)])
% ylim([-1.548,-1.536])
% ylim([-1.17196,-0.3343])
hold off;

%=========================== FIGURE SETTINGS =============================%
%% Figure resizing
ax_h = -0.08; ax_bot = 0; ann_pos = 0.45;    %<-------- Edit as necessary to fit figure properly
fig_width = 500; fig_height = 500;
ch = fig2.Children;

set(fig2, 'Position', [50, 50, fig_width, fig_height])
for par = 1:1:1
    
    % Plot points
    ax = ch(par);
    
    x = get(ax,'XTickLabel');
    set(ax,'XTickLabel',x,'FontName','Times','fontsize',14)
    y = get(ax,'XTickLabel');
    set(ax,'XTickLabel',y,'FontName','Times','fontsize',14)
    
    outerpos = ax.OuterPosition;

    ti = ax.TightInset;
    ti(1:2:3) = ti(1:2:3) - 0.015;
    ti(2:2:4) = ti(2:2:4) + 0.058;

    ax_width = outerpos(3) - ti(1) - ti(3) - 0.06;
    ax_height = outerpos(4) - ti(2) - ti(4) + ax_h;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2) - ax_bot;

    pos = [left bottom ax_width ax_height];
    ax.Position = pos;
    
end

legend_items = [est_h ds_h true_h];
legend_entries = [{'Estimated Position'}, {'Downsampled_position'}, {'True Position'}];

lh = legend(legend_items,legend_entries,'Orientation','vertical');
% lh = legend([h_obj lc_h],[title_obj label_lc],'Orientation','horizontal');
rect = [0.2613, 0.25, .27, .0528];
% rect = [0.3813, 0.85, .27, .0528];
set(lh, 'Position', rect, 'interpreter', 'latex', 'fontsize', 16)

if export_fig
    set(fig2,'color','w');
    savefig(fig2,'./Sample_plots/fig1.fig')
    export_fig(['./Sample_plots/1D_response_vel_delay_',num2str(n_delay),'_c_',num2str(c),'_l_',num2str(lambda),'.pdf'],'-p0.002',fig2); 
    export_fig(['./Sample_plots/1D_response_vel_delay_',num2str(n_delay),'_c_',num2str(c),'_l_',num2str(lambda),'.png'],'-p0.002','-r600',fig2); 
end

%% Plot algorithm output velocity

fig3 = figure('Name', 'Velocity');
hold on;
est_h = plot(time_raw, vel_est_cpp(:,2), '-r', 'linewidth',1.5);
% ds_h = plot(time_raw, groundtruth_pos_ds(:,2), '-g', 'linewidth',1.5);
true_h = plot(time_raw, groundtruth_vel_true(:,2), '-k', 'linewidth',1.5);
xlabel('Time (s)','interpreter','latex');
ylabel('y (m)','interpreter','latex');
xlim([min(time_raw),max(time_raw)])
% ylim([-1.548,-1.536])
% ylim([-1.17196,-0.3343])
hold off;

%=========================== FIGURE SETTINGS =============================%
%% Figure resizing
ax_h = -0.08; ax_bot = 0; ann_pos = 0.45;    %<-------- Edit as necessary to fit figure properly
fig_width = 500; fig_height = 500;
ch = fig3.Children;

set(fig3, 'Position', [50, 50, fig_width, fig_height])
for par = 1:1:1
    
    % Plot points
    ax = ch(par);
    
    x = get(ax,'XTickLabel');
    set(ax,'XTickLabel',x,'FontName','Times','fontsize',14)
    y = get(ax,'XTickLabel');
    set(ax,'XTickLabel',y,'FontName','Times','fontsize',14)
    
    outerpos = ax.OuterPosition;

    ti = ax.TightInset;
    ti(1:2:3) = ti(1:2:3) - 0.015;
    ti(2:2:4) = ti(2:2:4) + 0.058;

    ax_width = outerpos(3) - ti(1) - ti(3) - 0.06;
    ax_height = outerpos(4) - ti(2) - ti(4) + ax_h;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2) - ax_bot;

    pos = [left bottom ax_width ax_height];
    ax.Position = pos;
    
end

legend_items = [est_h true_h];
legend_entries = [{'Estimated Velocity'}, {'True Velocity'}];

lh = legend(legend_items,legend_entries,'Orientation','vertical');
% lh = legend([h_obj lc_h],[title_obj label_lc],'Orientation','horizontal');
rect = [0.2613, 0.25, .27, .0528];
% rect = [0.3813, 0.85, .27, .0528];
set(lh, 'Position', rect, 'interpreter', 'latex', 'fontsize', 16)

if export_fig
    set(fig3,'color','w');
    savefig(fig3,'./Sample_plots/fig1.fig')
    export_fig(['./Sample_plots/1D_response_vel_delay_',num2str(n_delay),'_c_',num2str(c),'_l_',num2str(lambda),'.pdf'],'-p0.002',fig3); 
    export_fig(['./Sample_plots/1D_response_vel_delay_',num2str(n_delay),'_c_',num2str(c),'_l_',num2str(lambda),'.png'],'-p0.002','-r600',fig3); 
end

%% Plot algorithm output acceleration

fig4 = figure('Name', 'Acceleration');
hold on;
est_h = plot(time_raw, acc_i_cpp(:,2), '-r', 'linewidth',1.5);
% ds_h = plot(time_raw, groundtruth_pos_ds(:,2), '-g', 'linewidth',1.5);
true_h = plot(time_raw, groundtruth_acc_true(:,2), '-k', 'linewidth',1.5);
xlabel('Time (s)','interpreter','latex');
ylabel('y (m)','interpreter','latex');
xlim([min(time_raw),max(time_raw)])
% ylim([-1.548,-1.536])
% ylim([-1.17196,-0.3343])
hold off;

%=========================== FIGURE SETTINGS =============================%
%% Figure resizing
ax_h = -0.08; ax_bot = 0; ann_pos = 0.45;    %<-------- Edit as necessary to fit figure properly
fig_width = 500; fig_height = 500;
ch = fig4.Children;

set(fig4, 'Position', [50, 50, fig_width, fig_height])
for par = 1:1:1
    
    % Plot points
    ax = ch(par);
    
    x = get(ax,'XTickLabel');
    set(ax,'XTickLabel',x,'FontName','Times','fontsize',14)
    y = get(ax,'XTickLabel');
    set(ax,'XTickLabel',y,'FontName','Times','fontsize',14)
    
    outerpos = ax.OuterPosition;

    ti = ax.TightInset;
    ti(1:2:3) = ti(1:2:3) - 0.015;
    ti(2:2:4) = ti(2:2:4) + 0.058;

    ax_width = outerpos(3) - ti(1) - ti(3) - 0.06;
    ax_height = outerpos(4) - ti(2) - ti(4) + ax_h;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2) - ax_bot;

    pos = [left bottom ax_width ax_height];
    ax.Position = pos;
    
end

legend_items = [est_h true_h];
legend_entries = [{'Estimated acceleration'}, {'True acceleration'}];

lh = legend(legend_items,legend_entries,'Orientation','vertical');
% lh = legend([h_obj lc_h],[title_obj label_lc],'Orientation','horizontal');
rect = [0.2613, 0.25, .27, .0528];
% rect = [0.3813, 0.85, .27, .0528];
set(lh, 'Position', rect, 'interpreter', 'latex', 'fontsize', 16)

if export_fig
    set(fig4,'color','w');
    savefig(fig4,'./Sample_plots/fig1.fig')
    export_fig(['./Sample_plots/1D_response_acc_delay_',num2str(n_delay),'_c_',num2str(c),'_l_',num2str(lambda),'.pdf'],'-p0.002',fig4); 
    export_fig(['./Sample_plots/1D_response_acc_delay_',num2str(n_delay),'_c_',num2str(c),'_l_',num2str(lambda),'.png'],'-p0.002','-r600',fig4); 
end
