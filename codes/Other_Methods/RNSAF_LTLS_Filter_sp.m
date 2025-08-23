function [EW,w] = RNSAF_LTLS_Filter_sp(Paras,weights)

        if nargin >1
            w = weights;
        else
            w = Paras.coeffs;                
        end
        un = Paras.input_sig;
        dn = Paras.y;
        NLMS_flag = Paras.NLMS_flag;
        % 
        % un: input signal
        % dn: desire signal
        % Paras: other init parameters
        % initialization
       %   initialization
        [L,N] = size(Paras.analysis);           % L is filter banks length, N is number of subbands
        H = Paras.analysis;                            % Analysis
        F = Paras.synthesis;                            % Synthesis        
        sigma_in = Paras.sigma_in;
        sigma_out = Paras.sigma_out;        
        num_iterations = length(un);
        x = zeros(L,1);               % Tapped-delay line of input signal (Analysis FB)
        d = zeros(L,1);               % Tapped-delay line of desired response (Analysis FB)
        z = zeros(L,1);               % Tapped-delay line of error signal (Synthesis FB)        
%         w = Paras.coeffs;
        M = length(w);
        U = zeros(M,N);           % Tapped-delay lines of adaptive subfilters
        AdaptStart = Paras.AdaptStart;
        mu = Paras.mu;
%         w = zeros(Paras.L_unknown,1);
        theta = sigma_out^2 / sigma_in^2;       % 根据论文定义的参数 \theta_i
        w_bar = [w; -sqrt(theta)];              % 扩展权重向量 \bar{w}
        h = Paras.w_true; % 未知系统,稀疏表示系数
        en = zeros(1,num_iterations);
        X = zeros(L, N);                        % 更新buffer
        EW = zeros(1,num_iterations);
%         alpha = ones(1,N)*1e-4;     % For NLMS
        alpha = 0.1;
%         D = Paras.D;
        ComputeEML = 1;
        % adaption
        for n = 1:num_iterations
            % 当前输入与期望信号
            x = [un(n); x(1:end-1)];  % Fullband input vector for band partitioning
            d = [dn(n); d(1:end-1)];  % Fullband desired response vector for band partitioning
            X = [x,X(:,1:end-1)];
            if ComputeEML == 1
                EW(n) = norm(w - h)^2 / norm(h)^2;                % System error norm (normalized)
%                 u = [un(n); u(1:end-1)];
%                 UDerr(n) = (b-w)'*u;                       % Undisturbed error 
            end
            % 使用 CMFB 进行子带分解
            if (mod(n,N)==0)              % Tap-weight adaptation at
            % lowest sampling rate, 如果需要按照每一次buffer来计算的话需要这样写               
                U1 = (H'*X)';             % Partitioning u(n) 
                U2 = U(1:end-N,:);
                U = [U1', U2']';          % subband data matrix
                dD = H'*d;                % Row vector, each element is decimated desired signal
                                          % 子带误差 e_iD 和梯度更新
                eD = dD - U'*w;                        
                if n >= AdaptStart
                    if  ~NLMS_flag 
                    w_bar_squared = sum(w_bar.^2); % 使用实际权重向量 w 的范数
%                   x_tilede_squared = sum(x_tilde.^2);
                    numerator = -(ones(M,1) * eD' .* U * w_bar_squared + w * (eD.^2)');
                    denominator = (w_bar_squared^2 * (sum(U.^2,1))'+ alpha * w_bar_squared * eD.^2);
                    grad =  numerator./ (ones(M,1)*denominator');  % 梯度拆解成分子分母分别计算再加
                    w = w - mu*sum(grad,2);
                    Paras.iter = Paras.iter + 1;                    
                    w_bar = [w; -sqrt(theta)];
                    else
                          w = w + U*(eD./(sum(U.*U)+alpha)')*mu;   % NLMS
                    end
                end               
                z = F * eD + z;                                       
                en(n-N+1:n) = z(1:N); 
                z = [z(N+1:end); zeros(N,1)];
            end             
            % 计算 NMSD
            Paras.EW = EW;                
        end
end