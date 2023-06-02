function [k0] = semi_causal_decision_feedback_detecion(theta,c,rho,rho_star,kmax,y,n,threshold_alpha,threshold_beta)

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
% algorithm = 4 - the semi-causal, equipped with decision-feedback detection scheme,  (15a,15c),
%
%Inputs:
% theta = [\hat{a}_1,...,\hat{a}_n]^T - the estimates of autoregressive  coefficients
% c     = [\hat{c}_1,...,\hat{c}_n]^T - the estimates of interpolation  coefficients
% rho - the estimate of variance of white driving noise, 
% rho_star - the estimate of leave-one-out signal interpolation errors,
% kmax - maximum length of detection alarm
% y - the measured signal
% n - the order of autoregression
% threshold_alpha = \mu_{\alpha}^2
% threshold_beta  = \mu_{\beta}^2
%
%Output:
% k0 - denotes number of corrupted samples
%
% Date: 08.2016
% © All rights reserved.
%##########################################################################

%#################################### Initialization
t                = n;
counter_prediction     = 0;
counter_interpolation  = 0;

x                = y(t:-1:t-n+1);
Q                = zeros(n,n);
A                = zeros(n,n);
A(1,:)           = theta';
A(2:end,1:end-1) = eye(n-1);
C                = zeros(n,1);
C(1)             = 1;
%####################################


%#################################### Detection loop
     for k=1:kmax+2*n
         
        %############################ Prediction (15a,15c)
        
        x     = A*x;
        Q     = A*Q*A'+rho*(C*C');
        sigma = C'*Q*C;
        e     = y(t+k) - C'*x;
        alpha = e^2/sigma;
        
        if alpha > threshold_alpha && k <= kmax
            counter_prediction = 0;
        else
            L = Q*C/sigma;
            x = x + L*e;
            Q = Q - sigma*(L*L');
            
            counter_prediction = counter_prediction + 1;
        end
        
        %############################ Interpolation (17)-(21)
        
        ym = 0;
        for ii=1:n
            ym = ym + c(ii)*(y(t+k-ii)+y(t+k+ii));    
        end

        e_star  = y(t+k) - ym;             
        beta    = (e_star^2)/(rho_star);
        
        if (beta > threshold_beta && k <= kmax) || k==1
            counter_interpolation = 0;
        else
            counter_interpolation = counter_interpolation + 1;
        end
        
        %############################
        
        if (counter_prediction == n) || (counter_interpolation == n)
            k0 = k - n;
            break;
        end
        %############################
        
     end
%####################################     
  
     
     if k >= kmax+n
         k0 = kmax;
     end
     
%####################################
end