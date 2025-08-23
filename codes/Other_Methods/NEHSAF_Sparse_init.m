function Paras = NEHSAF_Sparse_init(M,mu,lambda,flags,Paras)


    % Arguments:
    % M             Length of corresponding fullband filter
    % mu            Step size
    % N             Number of subbands
    % D             Decimation factor, D=N (Critical decimation), D<N (oversampling)

    % 参数初始化
  if  isfield(Paras, 'num_iterations')
  else
      Paras.num_iterations = 2e5;
  end
            % 迭代次数
    Paras.L_unknown = M;               % 未知系统的滤波器长度
    Paras.alpha = 0.1;                 % smooth factor
    Paras.mu = mu;                     % 迭代步长
%     Paras.sigma_in = 0;                % 输入/输出背景噪声标准差
%     Paras.p = 0.01;                    % 冲击噪声的概率
%     Paras.sigma_impulsive = sqrt(100); % 冲击噪声的标准差
    Paras.iter = 0;
    Paras.coeffs = zeros(M,1); %
    Paras.lambda = lambda;    
%     Paras.noise_type = noise_type;
    Paras.var_lambda_flag = flags.var_lambda_flag;
    Paras.NEHSAF_flag = flags.NEHSAF_flag;
%     Paras.snr = 30;
    
%     % 生成未知系统的目标权重向量 h    
%     rng('default')                     % for reproduction
%     h = randn(Paras.L_unknown, 1);     % 高斯随机分布
%     h = h / norm(h);                   % 归一化滤波器系数
%     Paras.h = h;
    
    
    Paras.AdaptStart    = M;           % Running effect of the analysis and adaptive 
%                                        % filter, minimum L + M
end
