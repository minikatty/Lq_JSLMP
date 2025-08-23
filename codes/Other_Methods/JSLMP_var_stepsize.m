clc;clear;
% close all
mc_num = 1;
%% data generation and initialation
for mc = 1:mc_num 
    N = 512;% length of sparse presentation
    L = 300; % length of received signal
    alpha = 1.5;gamma =2e-3; % impulsive noise level 4e-3
    % b = 1.2;
    % n = 11;
    iter_max = 4e5; % single run iterations, total iter = iter_max * length(k_Change).
    rng(3)
    H = sqrt(1/L)*randn(L,N);% measurement matrix
        % H=orth(H);
    K =30; % sparsity: nonzeros element
    sp = K/L; % sparsity degrees;
    x = full(sprandn(N,1,sp)); % sparse presentation    
    x_true = x./norm(x); % normalize 
    n = starnd(alpha,gamma,0,0,L); %impulsive disturbance
    s = H*x_true; % corrupted signal    
    y.y = s + n';% received signal
    y.iter_max = iter_max;
    %% signal config
    q =1;q1 = q;q2 = q;p = 1.8;
    eta1 = 0.006;%0.15; % regular parameters range 0:0.005:0.06;% 6
    eta2 = 0.006;%0.15; % varying with p, a bigger p corresponds to a smaller eat % 6
     %% computing the weights
    [MSD_JSLMP(mc,:),~] = JSLMP3(y,q1,q2,p,H,x_true,eta1,eta2,'normal');
end
    %% plot 
    figure
    semilogy(mean(MSD_JSLMP,1));