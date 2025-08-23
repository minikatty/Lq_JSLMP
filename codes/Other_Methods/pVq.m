clc;clear;close all
delete(gcp('nocreate'));
parpool(8);
%% data generation and initialation
N = 512;% length of sparse presentation
L = 300; 
alpha = 1.5;gamma = 4e-3; % impulsive noise level
k = 30;% nz component <N
sp = k/N; % sparsity degree
a = 0.11; 
b = 1.2;
n =11;
eta1 = logspace(log10(a), log10(b), n); % regular parameters range 0:0.005:0.06;
eta2 = eta1; % varying with p, a bigger p corresponds to a smaller eta
% p = 2;% length = 11 p = 1:0.1:2;
q1 = 0:0.1:1;
p = 1;
Len_q1 = length(q1);
mc = 48;%  MC counts %
% q2 = q1;
MSD_min = zeros(length(q1),mc); % temp msd matrix
% MSD_min_mc = MSD_min;
opt_eta1 = zeros(length(q1), mc);
opt_eta2 = zeros(length(q1), mc);
rng(1)

%% Monte Carlo 

parfor MC = 1:mc
    
    H = sqrt(1/L)*randn(L,N);% measurement matrix
    % H=orth(H);
    x = full(sprandn(N,1,sp)); % sparse presentation
    
    x = x./norm(x); % normalize
    
    n = starnd(alpha,gamma,0,0,L); %impulsive disturbance
    
    s = H*x; % corrupted signal
    
    y = s + n';% received signal
    
    % figure; plot(y,'b-'); plot the signal        
        
        q1 = 0:0.1:1;% length = 11
%         q2 = 0:0.1:1;
        
        for j = 1:Len_q1 % length = 11
            
%             q2 = q1(j); %[0,0.1,0.2,0.3,0.4]
            
            [opt_eta1(j,MC),opt_eta2(j,MC),~,MSD_min(j,MC)] = bestreg(y, q1(j), q1(j), p, H, x, eta1, eta2); %reg1,reg2, MSD,minimum
            
        end 
        
        disp(['[' datestr(now,13) '] This is ' int2str(MC) ' times.'])
%     MSD_min_mc = MSD_min_mc + MSD_min./mc;
%     opteta1 = opteta1+opt.eta1(j)/mc;
%     opteta2 = opteta2+opt.eta2(j)/mc;
    
end
save(['./results/pvq' num2str(10*p) '.mat'],'opt_eta1','opt_eta2','MSD_min');
% min_msd = min(min(MSD_min_mc));
% quit;
% min_index = find(MSD_min_mc= min_msd);

%% plot the best p&q
% [P,Q1]=meshgrid(p,q1);
% figure
% mesh(P, Q1,MSD_min_mc);
% figure;bar3(MSD_min_mc);
% lgMSD=10*log10(MSD_min_mc);
% figure;bar3(lgMSD);
% figure; imagesc(lambda1, lambda2,lgMSD);
% figure;contourf(lambda1,lambda2,lgMSD');