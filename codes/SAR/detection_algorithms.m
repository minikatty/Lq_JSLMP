function [y,d] = detection_algorithms(y,M,n,threshold_alpha,threshold_beta,kmax,algorithm)

%##########################################################################
%Author:      Marcin Cio砮k, PhD
%
%Affiliation: Gda駍k University of Technology,  
%             Faculty of Electronics, Telecommunications and Computer Science, Department of Automatic Control
%             street Narutowicza 11/12, 80-233 Poland
%
%E-mail:      marcin.ciolek@pg.gda.pl
%
%Conference:  ICASSP 2017
%Paper:       DETECTION OF IMPULSIVE DISTURBANCES IN ARCHIVE AUDIO SIGNALS
%
%ADAPTIVE DETECTION:
%
% 4 unidirectional detection algorithms are considered: 
% algorithm = 1 - the causal, equipped with open-loop detection scheme,
% algorithm = 2 - the causal, equipped with decision-feedback detection scheme,
% algorithm = 3 - the semi-causal, equipped with open loop detection scheme, 
% algorithm = 4 - the semi-causal, equipped with decision-feedback detection scheme.
%
%Inputs:
% y - the measured signal
% n - the order of autoregression
% M - segment length, where M = 2k + 1
% threshold_alpha = \mu_{\alpha}^2
% threshold_beta  = \mu_{\beta}^2
% kmax - maximum length of detection alarm
%
%Outputs:
% y - the processed signal
% d - the estimate of pulse location function
%
% Date: 08.2016
% ? All rights reserved.
%##########################################################################




%##########################################################################
if algorithm > 2
  semi_causal = true;
else
  semi_causal = false;  
end



N = length(y);


alarm_triggered    = false;
approved_samples   = 0;
alarms_counter     = 0;
k0                 = 0;

d  = zeros(N,1);
p  = zeros(n+1,1);
%##########################################################################


%##########################################################################
%##################################### the cosinusoidal window (10)

w = zeros(M,1); %窗函数
k = (M-1)/2; %

for l=-k:k
    w(l+k+1) = cos(pi*l/(2*(k+1))); %定义一个余弦窗，个人理解是：1.如果存在频域处理，则为了防止栅栏效应和频谱泄露；
                                                 %                                          2.仅仅是提供插值，是的数据平滑。
end

%  the normalizing constant (9), this is an exact value for the effective
%  width (the effective number of observations)
%  of the cosinusoidal window
L = k+1;                   %窗函数的归一化常数

%#####################################

%#####################################
% Results of following calculations are used next in the recursion of (8)
% %为了递归计算互相关矩阵的一些常数值存放
exp_value2_re = zeros(n+1,1);
exp_value2_im = exp_value2_re;
exp_value4_re = exp_value2_re;
exp_value4_im = exp_value2_re;
cos_value     = exp_value2_re;

exp_value1_re = real(exp(-1j*pi/(k+1)));
exp_value1_im = imag(exp(-1j*pi/(k+1)));
exp_value3_re = real(exp(1j*pi*k/(k+1)));
exp_value3_im = imag(exp(1j*pi*k/(k+1)));

for l=0:n

    exp_value2_re(l+1) = real(exp(1j*pi*l/(k+1)));
    exp_value2_im(l+1) = imag(exp(1j*pi*l/(k+1)));
    exp_value4_re(l+1) = real(0.5*exp(-1j*pi*l/(2*k+2)));   
    exp_value4_im(l+1) = imag(0.5*exp(-1j*pi*l/(2*k+2)));   
    cos_value(l+1)     = 0.5*cos(pi*l/(2*k+2));
    
end
%#####################################

%#####################################   
% Calculation of the quantities f_i(t) and g_i(t), which are recursively
% computable 

        t = M;
        
        f = zeros(n+1,1);
        g = zeros(n+1,1);
        h = zeros(N,n+1);
        
        for i=0:n
            
            ii = i+1;
            
            for l=-k+i:k
                
                h(t-k+l,ii) = y(t-k+l)*y(t-k+l-i);                             
                f(ii)       = f(ii) + h(t-k+l,ii);
                g(ii)       = g(ii) + h(t-k+l,ii)*exp(1j*pi*l/(k+1));
                  
            end
        end
       
       
        g_real = real(g);
        g_imag = imag(g);
        
        for i=0:n
            ii = i+1;           
                p(ii) = cos_value(ii)*f(ii) + g_real(ii)*exp_value4_re(ii) - g_imag(ii)*exp_value4_im(ii);  
        end 

%The  Yule-Walker (YW) estimates (6) can be obtained using the Levinson-Durbin 
%algorithm which guarantees model stability. 
        r = p/L;
        [theta,rho] =  Levinson_Durbin_algorithm(r,n);
%#####################################        
        
             
      
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%Main loop starts here


t_end = N-kmax-2*n;


for t=M+1:t_end
    

   
            if approved_samples <=0
                
            %##################################### the one-step-ahead prediction error (12)
            %#####################################  
            
            phi    = y(t-1:-1:t-n);                     
            e      = y(t) - phi'*theta;                                  
            alpha  = (e^2)/rho;   
            
            %#####################################
            %#####################################
           

            %##################################### prediction based detector (11)
            %#####################################
             if alpha > threshold_alpha %基本的冲击检测条件，3sigma准则变体，记作规则(1)
               
                if semi_causal == true %半因果增强机制
                     
                           
                %########################################################## interpolation based detector
                %##########################################################      
                % the leave-one-out signal interpolation (3)-(5)
                [c,rho_star] = leave_one_out_signal_interpolation(theta,rho,n);
                     
                %the strengthened alarm triggering rule (20)
                for tt=t-n+1:t

                    y_tilde = 0;
                    for i=1:n
                        y_tilde = y_tilde + c(i)*(y(tt-i)+y(tt+i));    %论文[18]中的一个特例，也就是t时刻缺损的值由前n个点和
                        % 后n个未受干扰的值来进行插值，观测信号其实就是[y(t-n)... y(t-1) y(t)
                        % y(t+1)...y(t+n)] 这里并不是进行的重构计算，这里只是增强检测算法必要的步骤。
                    end

                    e_star  = y(tt) - y_tilde;             
                    beta    = (e_star^2)/(rho_star);

                    if beta > threshold_beta     
                        alarm_triggered = true;
                        break;
                    end 
                
                end %增强检测机制到此停止，且增强检测是一种半因果的检测算法。它需要当前时刻之前与之后n个时刻的值
                
                %##########################################################
                %##########################################################
                
                 else
                     
                     alarm_triggered = true; %规则(1)下的警告触发
                     
                end %半因果判断检测结束
                 
             end    %基本准则(1)结束        
             
            end
            %#####################################
            %#####################################
        
            %##################################### detection and reconstruction schemes
            %#####################################
             if alarm_triggered == true
                 
          
                    
                     
                    %#################################################################
                    %#################################################################
                    yc = y(t-n:t+kmax+2*n);% 
                    
                    switch(algorithm) % 4 unidirectional detection algorithms %检测算法开启
                     
                         case 1
                                k0 = causal_open_loop_detecion(theta,rho,kmax,yc,n,threshold_alpha);
                                
                         case 2
                                k0 = causal_decision_feedback_detecion(theta,rho,kmax,yc,n,threshold_alpha);

                         case 3        
                                k0 = semi_causal_open_loop_detecion(theta,c,rho,rho_star,kmax,yc,n,threshold_alpha,threshold_beta);

                         case 4
                                k0 = semi_causal_decision_feedback_detecion(theta,c,rho,rho_star,kmax,yc,n,threshold_alpha,threshold_beta);

                    end
                     %k0 denotes number of corrupted samples
                    %#################################################################
                    %#################################################################               
      
                    % the AR-model based reconstruction of corrupted 
                    % samples can be carried out independently without any
                    % information loss for each local analysis frame starting and
                    % ending with n undistorted samples y(t)
                    t0 = t - 1;      
                    yo = vertcat(y(t0-n+1:t0), y(t0+k0+1:t0+k0+n));
                    y_tilde = AR_model_based_signal_interpolation(theta,n,k0,yo);                   
                    
                    y(t0+1:t0+k0)  = y_tilde;
                    d(t0+1:t0+k0)  = 1;
                    
                    approved_samples = k0 + n;

                    alarms_counter  = alarms_counter + 1;
                    alarm_triggered = false;
                    
                    %#################################################################
                    %#################################################################
             end
            %#####################################
            %#####################################
            

             approved_samples = approved_samples - 1;
        
           
            %#####################################
            %##################################### 
            %recursion of p_i(t) requires only  12 real add-multiply
            %operations per time update for a selected value of i and does not depend on M  
                             
            for i=0:n
                ii = i+1;

                h(t,ii)    = y(t)*y(t-i);                   
                f(ii)      = f(ii) - h(t-M+i,ii) + h(t,ii);

                G_re       = exp_value1_re*g_real(ii) - exp_value1_im*g_imag(ii)  + h(t-M+i,ii)*exp_value2_re(ii) + h(t,ii)*exp_value3_re;
                G_im       = exp_value1_im*g_real(ii) + exp_value1_re*g_imag(ii)  + h(t-M+i,ii)*exp_value2_im(ii) + h(t,ii)*exp_value3_im;
                g_real(ii) = G_re;
                g_imag(ii) = G_im;

                p(ii)      = cos_value(ii)*f(ii)  + G_re*exp_value4_re(ii) - G_im*exp_value4_im(ii);

            end

            %#####################################
            %#####################################

             %the Levinson-Durbin algorithm
             r = p/L;
             [theta,rho] =  Levinson_Durbin_algorithm(r,n); %系统参数的重新估计
            
            %#####################################
            %#####################################
    
end




end



               
                            

