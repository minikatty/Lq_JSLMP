function [y] = reconstruction_algorithm(y,n,M,d)
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
% Reconstruction algorithm
%
%Inputs:
% y - the measured signal
% n - the order of autoregression
% M - segment length, where M = 2k + 1
% d - the estimate of pulse location function
%
%Output:
% y - the reconstructed signal
%
% Date: 08.2016
% ? All rights reserved.
%##########################################################################


%##################################### the cosinusoidal window
w = zeros(M,1);
k = (M-1)/2;

for l=-k:k
    w(l+k+1) = cos(pi*l/(2*(k+1))); 
end

%  the normalizing constant (9), this is an exact value for the effective
%  width (the effective number of observations)
%  of the cosinusoidal window
L = k+1;
%#####################################

d_diff    = diff(d); 
alarm_beg = find(d_diff > 0.1)+1;
alarm_end = find(d_diff <-0.1);

N         = length(alarm_end);% the number of impulsive disturbance
rec_end   = length(y);

p = zeros(n+1,1); %Yule-Walker equation
%##########################################################################
    for t=1:N


            t0 = alarm_beg(t) - 1;

            if (t0-M+1 > 0) && (alarm_end(t) + n <= rec_end)


            %#####################################
            %#####################################
                k0  = alarm_end(t) - alarm_beg(t) + 1; %the length of impulsive disturbance
                yc  = y(t0-M+1:t0);
                yw  = yc.*w;

                for i=0:n  
                    temp = 0;
                    for l=1+i:M
                        temp = temp+ yw(l)*yw(l-i);
                    end     
                    p(i+1,1) = temp;
                end     

                r = p/L;
                [theta,~] = Levinson_Durbin_algorithm(r,n);

                % SIGNAL RECONSTRUCTION
                % the AR-model based reconstruction of corrupted 
                % samples can be carried out independently without any
                % information loss for each local analysis frame starting and
                % ending with n undistorted samples y(t)
                yo        = [y(t0-n+1:t0);y(t0+k0+1:t0+k0+n)];
                y_tilde   = AR_model_based_signal_interpolation(theta,n,k0,yo);

                y(t0+1:t0+k0) = y_tilde;

            %#####################################
            %#####################################


            end

    end
%##########################################################################

end



               
                            

