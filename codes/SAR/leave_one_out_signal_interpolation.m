
function [c,rho_star] = leave_one_out_signal_interpolation(a,rho,n)

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
% The leave-one-out signal interpolation (3)-(5)
%
% Numeration (x) corresponds to the equations in the paper
%
%Inputs:
% a   = [\hat{a}_1,...,\hat{a}_n]^T   - the estimates of the so-called autoregressive coefficients
% rho - the estimate of variance of white driving noise, 
% n   - the order of autoregression
%
%Outputs:
% c   = [\hat{c}_1,...,\hat{c}_n]^T   - the estimates of the interpolation coefficients
% rho_star - the estimate of leave-one-out signal interpolation errors,
%
% Date: 08.2016
% © All rights reserved.
%##########################################################################
  
    d  = 1 + a'*a;
 
  
   
    c  = zeros(n,1);
    %#############################
    for i=1:n-1

        temp = a(i);
        
        for j=1:n-i               
           temp = temp - a(j)*a(i+j);  
        end
        c(i) = temp/d;

    end
    %#############################
    c(n)     = a(n)/d;
    
    rho_star = rho/d;
    
end