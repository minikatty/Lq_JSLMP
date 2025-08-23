function [EW,w] = EHCAF_Sparse_Filter(H,x,y,Paras,weights)
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
        lambda = Paras.lambda;
%         alpha = ones(1,N)*1e-4;     % For NLMS
        alpha = Paras.alpha;
%         D = Paras.D;
        ComputeEML = 1;
%         zeta = 0.01/L;
        % adaption
        for n = 1:num_iterations
            % 当前输入与期望信号   
                k = mod(n,L)+1; 
                dn(n) = y(k);
                x = H(k,:)'; 
            
            if ComputeEML == 1
                EW(n) = norm(w - h)^2 / norm(h)^2;                % System error norm (normalized)
            end 
            
            en(n) = dn(n)-w'*x;
            incre = x * exp(-abs(cosh(lambda * en(n))^alpha)) * sinh(lambda * en(n))...
             * abs(cosh(lambda * en(n)))^(alpha-1) * sign(cosh(lambda * en(n)));
            w = w + mu*incre;
            
        end
 % 计算 NMSD
 % Paras.EW = EW;   
end