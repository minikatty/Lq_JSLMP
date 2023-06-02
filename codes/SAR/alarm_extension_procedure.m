function [d_extended] = alarm_extension_procedure(d,n,delta)
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
% Alarm extension:
%
% M. NiedŸwiecki and M. Cio³ek, “Elimination of impulsive disturbances from archive audio signals using 
% bidirectional processing,” Trans. Audio Speech Lang. Process., vol. 21, pp. 1046–1059,2013.
%
% Unlike artificially generated noise pulses,
% real impulsive disturbances corrupting audio signals are rarely
% confined to isolated samples. Moreover, most of them have
% “soft” edges (the more so, the higher sampling rate) which
% stems from the typical geometry of local damages of the
% recording medium (e.g., groove damages). The straightforward
% consequence of this fact is that detection alarms are seldom
% triggered at the very beginning of noise pulses. This may lead
% to small but audible distortions of the reconstructed audio
% material. Although detection delays can be reduced, or even
% eliminated, by lowering the detection multiplier , i.e., by
% making the outlier detector more sensitive to “unpredictable”
% signal changes, the improvement comes at a price: low detection
% thresholds may dramatically increase the number and
% length of detection alarms, causing the overall degradation of
% the results. An alternative solution, which works pretty well
% in practice, is based on shifting back the beginning of each
% detection alarm (once determined) by a small fixed number of
% samples.
% 
%Inputs:
% d - the estimate of pulse location function
% n - the order of autoregression
% delta - the shifting parameter
%
%Outputs:
% d_extended - the estimate of pulse location function after extensions
% 
%
% Date: 08.2016
% © All rights reserved.
%##########################################################################

%##########################################################################

d_differences    = diff(d);
alarm_beginning  = find(d_differences == 1)+1;
alarm_end        = find(d_differences == -1);

N                = length(d);
d_extended       = zeros(N,1);

number_of_alarms = length(alarm_end);


t2_previous = 0;
%##########################################################################
for i=1:number_of_alarms
    
                
    
                %##########################################################
                t1 = alarm_beginning(i);
                t2 = alarm_end(i);
                
                
                if i-1 > 0
                    t2_previous = alarm_end(i-1);
                end
                
                t1_new = t1;
                for j = 1:delta                    
                        if (t2_previous < t1 - j - n || t2_previous == 0)                                                      
                            t1_new = t1 - j;                           
                        end
                end
                             
                d_extended(t1_new:t2) = 1;
                %##########################################################

end
%##########################################################################


end