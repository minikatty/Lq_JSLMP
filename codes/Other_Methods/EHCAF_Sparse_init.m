function Paras = EHCAF_Sparse_init(M,mu,lambda,Paras)

    % Arguments:
    % M             Length of corresponding fullband filter
    % mu            Step size
    % N             Number of subbands
    % D             Decimation factor, D=N (Critical decimation), D<N (oversampling)

    % 参数初始化
  if  isfield(Paras, 'num_iterations')
  else
      Paras.num_iterations = 1.2e5; % 迭代次数
  end
            
    Paras.L_unknown = M;               % 未知系统的滤波器长度
    Paras.alpha = 6;                 %  factor
    Paras.mu = mu;                     % 迭代步长
    Paras.iter = 0;
    Paras.coeffs = zeros(M,1); %
    Paras.lambda = lambda;    
%     Paras.var_lambda_flag = flags.var_lambda_flag;
%     Paras.NEHSAF_flag = flags.NEHSAF_flag;
    
end
