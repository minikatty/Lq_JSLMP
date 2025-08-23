function [EW,w] = NEHSAF_Sparse_Filter(H,x,y,Paras,weights)
        % x is true known system
        % received signal： y
        % weights: initial system

        if nargin >4
            w = weights;
        else
            w = Paras.coeffs;                
        end
        %   initialization 
        num_iterations = Paras.num_iterations;
        [L,~] = size(H);               % 未知系统的滤波器长度
        dn = zeros(num_iterations,1);               % Tapped-delay line of input signal       
        
%         AdaptStart = Paras.AdaptStart;
        mu = Paras.mu;
%       w = zeros(Paras.L_unknown,1);       
        h = x; % 未知系统
        en = zeros(1,num_iterations);
        EW = zeros(1,num_iterations);
        var_lambda_flag = Paras.var_lambda_flag;
        NEHSAF_flag = Paras.NEHSAF_flag;
        lambda = Paras.lambda;
%         alpha = ones(1,N)*1e-4;     % For NLMS
        alpha = 0.98;
%         D = Paras.D;
        ComputeEML = 0; %0
        lambda_min = 0.1;
        lambda_max = 3;
        Expectation = 0; % for lambda mechanism
        beta = 0.1;
%         zeta = 0.01/L;
        % adaption
        for n = 1:num_iterations
            % 当前输入与期望信号   
                k = mod(n,L)+1; 
                dn(n) = y(k);
                x = H(k,:)'; 
%             else
%                 x = [un(n); x(1:end-1)];  %  input vector
%                 dn(n) = dn(n);
            
            if ComputeEML == 1
                EW(n) = norm(w - h)^2 / norm(h)^2;                % System error norm (normalized)
            end                

                if var_lambda_flag
                    Expectation =(1-beta)*Expectation + beta*exp(-sinh(lambda*en(n).^2));
                    lambda = alpha * lambda +(1-alpha) * (1-Expectation);
                    if lambda >lambda_max
                        lambda = lambda_max;
                    elseif lambda <lambda_min
                        lambda = lambda_max;
                    end
                end
             if ~NEHSAF_flag
                en(n) = dn(n)-w'*x;
                sinh2lambda_en = sinh(2 * lambda * en(n));
                incre = mu * exp(-sinh2lambda_en^2) * sinh2lambda_en * x;   %(EHSAF)
                w = w +incre;

             else %(NEHSAF)
                en(n) = dn(n)-w'*x;
                norm_x = norm(x);
                normal_grad =  exp(-sinh(lambda * en(n)/norm_x)^2) * sinh(2*lambda*en(n)/norm_x) * x/norm_x;     % + zeta*sign(x)        
                w = w + mu*normal_grad;                 
             end                
        end
 % 计算 NMSD
 % Paras.EW = EW;   
end