function  [y_f,y_fb,d_fb] = bidirectional_signal_processing(y,n,M,threshold_alpha,threshold_beta,kmax,algorithm)
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
% 4 unidirectional detection algorithms are considered: 
% algorithm = 1 - the causal, equipped with open-loop detection scheme,
% algorithm = 2 - the causal, equipped with decision-feedback detection scheme,
% algorithm = 3 - the semi-causal, equipped with open loop detection scheme, 
% algorithm = 4 - the semi-causal, equipped with decision-feedback detection scheme,
% and 4 bidirectional algorithms i.e.,  noncausal extensions of above algorithms.
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
% y_f  - the restored signal in forward-time
% y_fb - the restored signal, bidirectional signal processing
%
% Date: 08.2016
% ? All rights reserved.
%##########################################################################

%##########################################################################
%##########################################################################

display('forward-time signal processing')
[y_f,d_f] = detection_algorithms(y,M,n,threshold_alpha,threshold_beta,kmax,algorithm);

% Noncausal detection
% M. NiedŸwiecki and M. Cio³ek, “Elimination of impulsive disturbances from archive audio signals using 
% bidirectional processing,? Trans. Audio Speech Lang. Process., vol. 21, pp. 1046?1059,2013.
% 
display('backward-time signal processing')
[y_b,d_b] = detection_algorithms(flipud(y),M,n,threshold_alpha,threshold_beta,kmax,algorithm);    

display('Alarm extension')
delta     = 2;
[d_ext_f] = alarm_extension_procedure(d_f,n,delta);
[d_ext_b] = alarm_extension_procedure(d_b,n,delta);
          

display('Combine forward-time and backward-time detection_alarms based on the set of local fusion rules')
[d_local_fusion_rules,~,~,~] =  combine_forward_backward_detection_alarms(d_ext_f,flipud(d_ext_b),d_f,flipud(d_b),n,delta,kmax);

d_fb = d_local_fusion_rules;
         
display('forward-time signal reconstruction')
[y_fb] = reconstruction_algorithm(y,n,M,d_fb);

%##########################################################################
%##########################################################################



end