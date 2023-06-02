function [k0] = causal_open_loop_detecion(theta,rho,kmax,y,n,threshold_alpha)

%##########################################################################
%Author:      Marcin Cio³ek, PhD
%
%Affiliation: Gdañsk University of Technology,  
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
% algorithm = 1 - the causal, equipped with open-loop detection scheme (15a,15b),
%
%Inputs:
% theta = [\hat{a}_1,...,\hat{a}_n]^T - the estimates of autoregressive  coefficients
% rho - the estimate of variance of white driving noise, 
% kmax - maximum length of detection alarm
% y - the measured signal
% n - the order of autoregression
% threshold_alpha = \mu_{\alpha}^2
%
%Output:
% k0 - denotes number of corrupted samples
%
% Date: 08.2016
% © All rights reserved.
%##########################################################################

%#################################### Initialization
t                   = n;
counter_prediction  = 0;

x                = y(t:-1:t-n+1);
Q                = zeros(n,n);
A                = zeros(n,n);
A(1,:)           = theta';
A(2:end,1:end-1) = eye(n-1);
C                = zeros(n,1);
C(1)             = 1;
%####################################


%#################################### detection loop
     for k=1:kmax+n
         
        %############################ (15a,15b))
        
        x     = A*x;
        Q     = A*Q*A'+rho*(C*C');
        sigma = C'*Q*C;
        e     = y(t+k) - C'*x;
        alpha = e^2/sigma;
        
        %############################
        
        if alpha > threshold_alpha && k <= kmax
            counter_prediction = 0;
        else
            counter_prediction = counter_prediction + 1;
        end
        
        if counter_prediction == n
            k0 = k - n;
            break;
        end
        %############################
        
     end
%####################################     
  
     
     if k == kmax+n
         k0 = kmax;
     end
     
%####################################
end