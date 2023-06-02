function [theta,rho] =  Levinson_Durbin_algorithm(r,n)

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
%ESTIMATION OF SIGNAL PARAMETERS: 
%The  Yule-Walker (YW) estimates (6) can be obtained using the Levinson-Durbin 
%algorithm which guarantees model stability. 
%T. Söderströmand and P. Stoica, System Identification. Englewoods Cliffs NJ: Prentice-Hall, 1988.

% Numeration (x) corresponds to the equations in the paper
%
%Inputs:
% n - the order of autoregression
% r = [\hat{r}_0,..., \hat{r}_n] - the vector of the local estimates of autocorrelation coefficients of y(t) (7) 
%
%Outputs:
% theta and rho are the YW estimates (6), where theta = [\hat{a}_1,...,\hat{a}_n]
%
% Date: 08.2016
% © All rights reserved.
%##########################################################################

            rho       = r(1);

            theta     = zeros(n,1);
            theta_old = theta;
            
            %#########################################
            %#########################################
            for k=1:n
                
                i=k+1;
                temp = 0;
    
                for j=1:k-1
                    temp =  temp + theta_old(j)*r(i-j);
                end
    
                    refl_coeff    = (r(i)- temp)/rho;
                    theta(k)      =  refl_coeff;
   
                for j=1:k-1
                    theta(j) = theta_old(j) - refl_coeff*theta_old(k-j);
                end
                theta_old    = theta;
                              
                rho = (1-refl_coeff^2)*rho;
            end         
            %#########################################
            %#########################################
            
end