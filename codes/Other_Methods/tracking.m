function [convg_MSD] = tracking(signal)
    algorithms = signal.algorithms;
    k_Change = signal.k_Change;
    num_k_change = length(k_Change);
    M =length(signal.true)/ num_k_change;    
    L = size(signal.H, 1);
    eta1 = signal.eta1;
    eta2 = signal.eta2;
    q1 = signal.q1;
    q2 = signal.q2;
    p = signal.p;   
    y.iter_max = signal.iter_max;
    convg_MSD = []; % zeros(1,num_k_change*signal.iter_max);
    
    for i = 1 : num_k_change
        
        H = signal.H(:,(i-1)*M+1:i*M);
         
        y.y = signal.input((i-1)*L+1:i*L);
         
        w_true = signal.true((i-1)*M+1:i*M);
         
        w_true =  w_true';
        
        switch algorithms
            case 'JSLMP'
                if i == 1
                    [MSD, tmp_weights] = JSLMP3(y,q1,q2,p,H,w_true,eta1,eta2,'normal');
                else 
                    [MSD, tmp_weights] = JSLMP3(y,q1,q2,p,H,w_true,eta1,eta2,'normal',tmp_weights);                      
                end
                convg_MSD = [convg_MSD, MSD'];
            case 'NEHSAF'
                % parameters config
                mu =  0.02;lambda = 0.5;
                Paras.num_iterations = y.iter_max;
                flags.var_lambda_flag = true;
                flags.NEHSAF_flag = true;
                Paras = NEHSAF_Sparse_init(M,mu,lambda,flags,Paras);
                if i == 1
                    [MSD, tmp_weights] = NEHSAF_Sparse_Filter(H,w_true,y.y,Paras); % NEHSAF_Sparse_Filter(H,x,y,Paras)
                else 
                    [MSD, tmp_weights] = NEHSAF_Sparse_Filter(H,w_true,y.y,Paras,tmp_weights);                      
                end
                convg_MSD = [convg_MSD, MSD];
            case 'RNSAF_LTLS'
                % parameters config
                N = 4;                       % Number of subbands
                mu = 0.2;                   % stepsize
%                 M = 512;                     % Length of corresponding fullband filter
                D = N/2;                     % Decimation factor, D=N (Critical decimation), D<N (oversampling)
                K = 4;                       % overlapping factor K=4
                L = 2*K*N;                   % Length of analysis filters, L = 2KN
                Paras.num_iterations = y.iter_max;
                Paras.NLMS_flag   = 1;
                idx = 1+(i-1)*y.iter_max:i*y.iter_max;
                Paras.n_RNSAF = signal.n_RNSAF(idx);
                Paras = RNSAF_LTLS_sp_init(H,w_true,M,mu,N,D,L,Paras);  % input is changed, little difference
            
                if i == 1
                     [MSD, tmp_weights] = RNSAF_LTLS_Filter_sp(Paras); % NEHSAF_Sparse_Filter(H,x,y,Paras)
                else 
                     [MSD, tmp_weights] = RNSAF_LTLS_Filter_sp(Paras,tmp_weights);                      
                end
                convg_MSD = [convg_MSD, MSD];
            case 'EHCAF'
                Paras.num_iterations = y.iter_max;
                lambda = 1.1;
                mu = 0.1;                
                Paras = EHCAF_Sparse_init(M,mu,lambda,Paras);
                if i == 1
                    [MSD, tmp_weights] = EHCAF_Sparse_Filter(H,w_true,y.y,Paras);
                else
                    [MSD, tmp_weights] = EHCAF_Sparse_Filter(H,w_true,y.y,Paras,tmp_weights);
                end
                convg_MSD = [convg_MSD, MSD];
        end
    end
end