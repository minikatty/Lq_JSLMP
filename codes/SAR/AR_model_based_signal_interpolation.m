function ym = AR_model_based_signal_interpolation(theta,n,k0,yo)
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
% AR model based signal interpolation:
% M.Nied焪iecki,揝tatistical reconstruction of multivariate time series,� 
% IEEE Trans. Signal Process., vol. 41,pp. 451�457, 1993.
% 
%
%Inputs:
% theta = [\hat{a}_1,...,\hat{a}_n]^T - the estimates of autoregressive  coefficients
% n - the order of autoregression
% k0 denotes number of corrupted samples
% yo = [y(t0-n+1:t0); y(t0+k0+1:t0+k0+n)]^T denotes vector of n samples preceding  and 
% n samples succeding the  block of corrupted samples
%
%Outputs:
% ym =[\hat{y}(t0+1),...,\hat{y}(t0+k0)]^T - the AR model based signal interpolation
%
% Date: 08.2016
% � All rights reserved.
%##########################################################################

    a = [1; -theta];
    A = zeros(n+k0,2*n+k0);

    for j=1:n+k0
        for i=1:n+1
            A(j,i+j-1) = a(n+2-i); 
        end
    end  
  
    Am = A(:,n+1:n+k0);
    Ao = horzcat(A(:,1:n), A(:,n+k0+1:end));

    ym = -(Am'*Am)\(Am'*Ao*yo);%-pinv(Am'*Am)*(Am'*Ao*yo);
    
end