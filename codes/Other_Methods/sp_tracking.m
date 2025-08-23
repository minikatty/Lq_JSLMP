clc;clear;close all
mc_num = 1;
%% data generation and initialation
for mc = 1:mc_num 
    N = 512;% length of sparse presentation
    L = 300; 
    alpha = 1.5;gamma =4e-3; % impulsive noise level 4e-3
    k_Change = [30]; % the number of sparse signal change k = 30;% nz component <N; 30 20 40
    % b = 1.2;
    % n = 11;
    eta1 = 0.15; % regular parameters range 0:0.005:0.06;% 6
    eta2 = 0.15; % varying with p, a bigger p corresponds to a smaller eat % 6
    iter_max = 1.5e5; % single run iterations, total iter = iter_max * length(k_Change).
    % MSD_min = zeros(1,length(q1)); % temp msd matrix
    % MSD_min_mc = MSD_min;
    % opteta1 = zeros(length(q1),1);
    % opteta2 = opteta1;opt_eta1=opteta1;opt_eta2 = opteta1;
    % rng(1)
    input_final = [];
    x_true = [];
    H_trans = [];
    n_RNSAF =[];

    %%  sythesis signal with changing sparsity
    for  i = 1: length(k_Change)   
        H = sqrt(1/L)*randn(L,N);% measurement matrix
        % H=orth(H);
        H_trans = [H_trans, H];
        sp = k_Change(i)/N; % sparsity degrees
        x = full(sprandn(N,1,sp)); % sparse presentation    
        x = x./norm(x); % normalize 
        rng(i + mc);
        n = starnd(alpha,gamma,0,0,L); %impulsive disturbance
        noise_level = [1.8,1.5,2.2];
        n_RNSAF2 = noise_level(i)*starnd(alpha,gamma,0,0,iter_max);
        n_RNSAF = [n_RNSAF n_RNSAF2];
        s = H*x; % corrupted signal    
        y = s + n';% received signal    
        x_true = [x_true x'];
        input_final = [input_final y'];
    end

    %% signal config
        q = 0.4;
        q1 = q;
        q2 = q;
        p = 1.1;
        sig.input = input_final;
        sig.iter_max = iter_max; % predefined iterations maximum every change
        sig.q1 = q1;
        sig.q2 = q2;
        sig.p = p;
        sig.H = H_trans; %
        sig.k_Change = k_Change;
        sig.true = x_true; 
        sig.eta1 = eta1;
        sig.eta2 = eta2;
        sig.n_RNSAF = n_RNSAF;
        % snr= 20*log10(sum(abs(s))/sum(abs(n')));
     %% computing the weights
        algorithms = 'JSLMP'; % include 4 algorithms: JSLMP, NEHSAF RNSAF_LTLS，and EHCAF.
        sig.algorithms = algorithms;
        MSD_JSLMP (mc,:)= tracking(sig);
        % sig.algorithms = 'NEHSAF';
        % MSD_NEHSAF(mc,:) = tracking(sig); 
        % sig.algorithms = 'RNSAF_LTLS';
        % MSD_RNSAF_LTLS(mc,:) = tracking(sig);
        % sig.algorithms = 'EHCAF';
        % MSD_EHCAF(mc,:) = tracking(sig);

end
    %% plot 
    figure
    semilogy(mean(MSD_JSLMP,1)); hold on;
    % semilogy(mean(MSD_NEHSAF,1)); hold on;
    % semilogy(mean(MSD_RNSAF_LTLS,1)); hold on;
    % semilogy(mean(MSD_EHCAF,1));
% for i = 1:length(k_Change)
%     10*log10(MSD(1+(i-1)*iter_max))
% end