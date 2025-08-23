function Paras = RNSAF_LTLS_sp_init(H_trans,w_true,M,mu,N,D,L,Paras)    

    input_sig = reshape(H_trans', [], 1);   

    Paras.w_true = w_true;
    Paras.sigma_in = sqrt(0.0001);       % 输入/输出背景噪声标准差
    Paras.sigma_out = sqrt(0.0001);
    
        
%     [quotient, remainder] = idivide(num_iterations, length(output_sig));

    % Arguments:
    % M             Length of corresponding fullband filter
    % mu            Step size
    % N             Number of subbands
    % D             Decimation factor, D=N (Critical decimation), D<N (oversampling)
    % L             Length of analysis filters, M=2KN, K is overlapped factor
    % NLMS_flag = 0;
    % 参数初始化
    if  isfield(Paras, 'num_iterations')
    else
      Paras.num_iterations = 2e5;      % 迭代次数
    end   
    
    num_iterations = Paras.num_iterations;
    len_input = num_iterations+1-length(w_true);
    if length(input_sig) > len_input
        input_sig = input_sig(1:len_input); 
    end
    quotient = floor(len_input/length(input_sig));
    remainder = mod(len_input,length(input_sig));
%     tmp = conv(input_sig',w_true');
    % 循环使用数据
           
   if quotient > 2
        for i = 1:(quotient-1)
            input_sig = [input_sig;input_sig];
        end
   end
   if remainder~=0
        input_sig = [input_sig;input_sig(1:remainder)]; 
    %         output_sig = [output_sig;output_sig(1:remainder)];
   end

    output_sig = (conv(input_sig',w_true') + Paras.n_RNSAF)';
%     truncate =  ;
    Paras.input_sig = input_sig;
    Paras.y = output_sig;
%     Paras.L_unknown = M;               % 未知系统的滤波器长度
    Paras.alpha = 0.1;                 % 标度参数
    Paras.mu = mu;                     % 迭代步长
%     Paras.sigma_in = sqrt(0.01);       % 输入/输出背景噪声标准差
%     Paras.sigma_out = sqrt(0.01);
%     Paras.p = 0.01;                    % 冲击噪声的概率
%     Paras.sigma_impulsive = sqrt(100); % 冲击噪声的标准差
    Paras.N = N;                       % Number of subbands
    Paras.L = L;                       % Length of analysis filters, M=2KN,
    Paras.D = D;
    Paras.iter = 0;
    Paras.coeffs = zeros(M,1); % only N/2+1 subfilters are considered if DFT Filter banks used
    
%     % 生成未知系统的目标权重向量 h    
%     rng('default')                     % for reproduction
%     h = rand(Paras.L_unknown, 1) - 0.5;% 均匀分布 [-0.5, 0.5]
%     h = h / norm(h);                   % 归一化滤波器系数
%     Paras.h = h;
%     signal.H = signal.H;

%   tmp = reshape(H',[],1); % 重新生成desired数据
    % 生成分析、合成滤波器cosine-modulated filter banks
    L = Paras.L;
    N = Paras.N; 
    hopt = opt_filter(L-1,N);          % Create prototype lowpass filter
    [H,F] = make_bank(hopt,N);         % Generate  H and  F  cosine-modulated filter bank
    H = sqrt(N)*H';  F = sqrt(N)*F';   % Scaling
    Paras.analysis = H;                % Analysis
    Paras.synthesis = F;               % Synthesis
%     Paras.NLMS_flag = NLMS_flag;
    Paras.AdaptStart    = L + M;           % Running effect of the analysis and adaptive 
                                       % filter, minimum L + M
end
