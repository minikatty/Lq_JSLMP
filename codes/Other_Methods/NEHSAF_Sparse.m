clc;clear;close all
mc = 1;
%% data generation and initialation
N = 512;% length of sparse presentation,i.e. unknown system length;
L = 512; 
alpha = 1.5;gamma = 5e-3; % impulsive noise level
k = 30;% nz component <N
sp = k/N; % sparsity degree
% rng('Default')
% M = 16; % unknown system length;
mu = 0.02;
lambda = 0.5;
flags.var_lambda_flag = false;
flags.NEHSAF_flag = true;
% noise_type = "Bernouli"; % 'S_alpha_S', 'Bernouli',"Gaussian","Lp-Lq";
%% Monte Carlo 

for MC = 1:mc
    
    H = sqrt(1/L)*randn(L,N);% measurement matrix
    % H=orth(H);
    x = full(sprandn(N,1,sp)); % sparse presentation
    
    x = x./norm(x); % normalize
    
    n = starnd(alpha,gamma,0,0,L); %impulsive disturbance
    
    s = H*x; % clean signal
    
    y = s + n';% received signal
    
    Paras = NEHSAF_Sparse_init(N,mu,lambda,flags);
    
    [en,result] = NEHSAF_Sparse_Filter(H,x,y,Paras);
       
    err_sqr = en.^2;
    EW = result.EW;
end

%% 绘制收敛曲线
iter =  result.num_iterations;
figure;                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
% MSE = err_sqr;
hold on; plot((0:length(MSE)-1),10*log10(MSE));%/1024
axis([0 iter -200 50]); %/1024
xlabel('Number of iterations'); % (\times 1024 input samples)
ylabel('Mean-square error (with delay)'); grid on;

figure;                          % Plot misalignment
% q = 0.99; EW = filter((1-q),[1 -q],EW);
hold on; plot((0:length(EW)-1),10*log10(EW)); % /1024
xlabel('Number of iterations'); % (\times 1024 input samples) 
ylabel('Misalignment (dB)');
grid on;