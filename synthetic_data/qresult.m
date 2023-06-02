%% figure;
clear all;
clc;
load('qres.mat','qres');
iter_max1=6e4;
qres1=qres(:,1:iter_max1);
 iteration=(1:1:iter_max1);
 figure
%  space=8000; %MarkerSpace
semilogy(iteration,qres1(2,:),'-','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,qres1(6,:),'-*','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,qres1(5,:),'-^','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,qres1(1,:),'-o','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,qres1(4,:),'-d','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,qres1(3,:),'-s','LineWidth',1.5,'markersize',6);  grid on;
% set(legend,'fontname','Arial');
% set(legend,'fontsize',12);
str={'q=0', 'q=0.2', 'q=0.4','q=0.6', 'q=0.8' ,'q=1'};
% columnlegend(3,str);
%%
clear all;
clc;
load('presult.mat','p_result');
iter_max1=size(p_result,2);
qres1=p_result(:,1:iter_max1);
 iteration=(1:1:iter_max1);
 figure
%  space=8000; %MarkerSpace
semilogy(iteration,p_result(3,:),'-','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,p_result(4,:),'-*','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,p_result(5,:),'-^','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,p_result(6,:),'-o','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,p_result(7,:),'-d','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,p_result(8,:),'-s','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,p_result(1,:),'--','LineWidth',1.5,'markersize',6); hold on;
semilogy(iteration,p_result(2,:),'--*','LineWidth',1.5,'markersize',6);  grid on;
% set(legend,'fontname','Arial');
% set(legend,'fontsize',12);
str={'p=1', 'p=1.2', 'p=1.4','p=1.6', 'p=1.8' ,'p=2','p=1 (\mu=0.001)','p=2 (\mu=0.1)'};
columnlegend(3,str);
%%












