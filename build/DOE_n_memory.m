clc
clear all
close all
format compact
format long
addpath('Support_functions')
addpath('Support_functions\export_fig')

n_postrack = 20; % downsampling frequency in manipulated signal
n_delay = 20; % to label the figure
% n_memory_vec = [0:1:n_postrack];
n_memory_vec = [0:1:61];
lambda = 2;
c = 15;

ll_c = 0.0; ul_c = 100;
ll_lambda = 0.0; ul_lambda = 50;
lb = [ll_c ll_lambda]; % lower bounds for lambda , C
ub = [ul_c ul_lambda]; % upper bounds for lambda , C

X = [c, lambda];
X = scaling(X,lb,ub,1); % scale variables between 0 and 1

for i=1:1:length(n_memory_vec)
        n_memory = n_memory_vec(i);
        
        param = {1,lb,ub,n_postrack,n_memory};
        [Err(i),v_Err(i),Fcont(i)] = RMSF_cpp(X',param);
end

%% Plot RMS results
fig1 = figure(1);
plot(n_memory_vec,Err,'-r','linewidth',2)
xlabel('$n_{memory}$','fontsize',14,'interpreter','latex')
ylabel('$RMS_{error}$','fontsize',14,'interpreter','latex')
% xlim([n_memory_vec(1),n_memory_vec(end)])

hold on
% mark minimum position
[err_min,n] = min(Err);
min_n_memory=n_memory_vec(n); %your point goes here 
line([min_n_memory min_n_memory],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','linewidth',1.5)

x = get(gca,'XTickLabel');
set(gca,'XTickLabel',x,'FontName','Times','fontsize',18)
set(gca,'XTickLabelMode','auto')
y = get(gca,'YTickLabel');
set(gca,'YTickLabel',y,'FontName','Times','fontsize',18)

ann_pos = 0.55;
A = [ann_pos 0.6 0.1 0.1];
t = annotation('textbox',A,'String',['anti-delay = ',num2str(min_n_memory),' units'],'LineStyle','none','fontsize',14);
set(t,'interpreter','latex')

set(fig1,'color','w');
export_fig(['./Sample_plots/DOE_n_memory_delay_',num2str(n_delay),'.pdf'],'-p0.002',fig1);
export_fig(['./Sample_plots/DOE_n_memory_delay_',num2str(n_delay),'.png'],'-p0.002','-r600',fig1); 

%% Plot velocity RMS results
fig2 = figure(2);
plot(n_memory_vec,v_Err,'-r','linewidth',2)
xlabel('$n_{memory}$','fontsize',14,'interpreter','latex')
ylabel('Velocity $RMS_{error}$','fontsize',14,'interpreter','latex')
% xlim([n_memory_vec(1),n_memory_vec(end)])

hold on
% mark minimum position
[err_min,n] = min(v_Err);
min_n_memory=n_memory_vec(n); %your point goes here 
line([min_n_memory min_n_memory],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','linewidth',1.5)

x = get(gca,'XTickLabel');
set(gca,'XTickLabel',x,'FontName','Times','fontsize',18)
set(gca,'XTickLabelMode','auto')
y = get(gca,'YTickLabel');
set(gca,'YTickLabel',y,'FontName','Times','fontsize',18)

ann_pos = 0.55;
A = [ann_pos 0.6 0.1 0.1];
t = annotation('textbox',A,'String',['anti-delay = ',num2str(min_n_memory),' units'],'LineStyle','none','fontsize',14);
set(t,'interpreter','latex')

set(fig2,'color','w');
export_fig(['./Sample_plots/DOE_n_memory_delay_VRMS_',num2str(n_delay),'.pdf'],'-p0.002',fig2);
export_fig(['./Sample_plots/DOE_n_memory_delay_VRMS_',num2str(n_delay),'.png'],'-p0.002','-r600',fig2); 

%% Plot FCONT results
fig3 = figure(3);
plot(n_memory_vec,Fcont,'-r','linewidth',2)
xlabel('$n_{memory}$','fontsize',14,'interpreter','latex')
ylabel('$|f|$','fontsize',14,'interpreter','latex')
% xlim([n_memory_vec(1),n_memory_vec(end)])

hold on
% mark minimum position
[err_min,n] = min(Fcont);
min_n_memory=n_memory_vec(n); %your point goes here 
line([min_n_memory min_n_memory],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','linewidth',1.5)

x = get(gca,'XTickLabel');
set(gca,'XTickLabel',x,'FontName','Times','fontsize',18)
set(gca,'XTickLabelMode','auto')
y = get(gca,'YTickLabel');
set(gca,'YTickLabel',y,'FontName','Times','fontsize',18)

ann_pos = 0.55;
A = [ann_pos 0.6 0.1 0.1];
t = annotation('textbox',A,'String',['anti-delay = ',num2str(min_n_memory),' units'],'LineStyle','none','fontsize',14);
set(t,'interpreter','latex')

set(fig3,'color','w');
export_fig(['./Sample_plots/DOE_n_memory_delay_FCONT',num2str(n_delay),'.pdf'],'-p0.002',fig3);
export_fig(['./Sample_plots/DOE_n_memory_delay_FCONT',num2str(n_delay),'.png'],'-p0.002','-r600',fig3); 
%% Save the DOE results
save('DOE_results_n_memory','n_memory_vec','Err','v_Err','Fcont')