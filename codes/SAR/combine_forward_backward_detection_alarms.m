function [d_local_fusion_rules,d_global_union_rules,d_global_intersecion_rules,pattern_counter] = ...
    combine_forward_backward_detection_alarms(d_f,d_b,d_f_org,d_b_org,n,delta,kmax)
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
% Bidirectional Detection of Noise Pulses:
%
% M. NiedŸwiecki and M. Cio³ek, “Elimination of impulsive disturbances from archive audio signals using 
% bidirectional processing,? Trans. Audio Speech Lang. Process., vol. 21, pp. 1046?1059,2013.
%
% The simplest approach to combining results of forward-time
% and backward-time detection is the one based on global decision
% rules, such as the intersection rule or the union rule.
% In the first case detection alarm is raised only when the sample
% is questioned by both detectors, and in the second case—when it
% is questioned by at least one of the detectors. Preliminary tests
% have shown that neither of these rules works satisfactorily in
% practice. The intersection rule is too conservative—it tends to
% overlook many small noise pulses and produces underfitted (too
% short) detection alarms. The union rule is too liberal—it yields
% many overfitted (too long) detection alarms which, after interpolation,
% result in audible signal distortions.
% To avoid problems mentioned above, different configurations
% of forward and backward detection alarms were divided into several classes and subclasses.
% Each class was analyzed separately in order to determine
% the best way of combining detection alarms. The final detection
% decision is a result of application of a certain number of
% local, case-dependent decision rules, called atomic fusion rules,
% rather than using a single global rule applicable to all cases.
% 
%Inputs:
% d_f     - the forward-time estimate of pulse location function after  extension
% d_b     - the backward-time estimate of pulse location function after  extension
% d_f_org - the forward-time estimate of pulse location function 
% d_b_org - the backward-time estimate of pulse location function 
% n       - the order of autoregression
% delta   - the shifting parameter
% kmax    - maximum length of detection alarm in unidirectional processing
%
%Outputs:
% d_local_fusion_rules       - combining results of forward-time  and backward-time detection is based on atomic fusion rules
% d_global_union_rules       - combining results of forward-time  and backward-time detection is based on global union rule
% d_global_intersecion_rules - combining results of forward-time  and backward-time detection is based on global intersection rule
%
% Date: 08.2016
% ? All rights reserved.
%##########################################################################

%#############################################################
%############################################################# merging
% we will sort out detection alarms in consecutive analysis
% frames defined as the minimum-length intervals that start and end with n no-alarm decisions

max_frame       = 2*kmax; %maximum length of detection alarm is doubled

d_fb            = d_f | d_b;

d_fb_diff       = diff(d_fb);
d_fb_beginning  = find(d_fb_diff >  0.1) +1;
d_fb_end        = find(d_fb_diff < -0.1);


analysis_frames = d_fb;

L = length(d_fb_end);

for i=2:L
   
    if d_fb_beginning(i) - d_fb_end(i-1)  <= n        
        analysis_frames(d_fb_end(i-1):d_fb_beginning(i)) = 1;
        % if the length between the next impulsive detection and previous
        % detection end is less than the orders of AR model, the data
        % within the range is all set to 1, i.e., regarded as outliners 
    end
    
end %(maximum no-alarm frames principle, alarm sepatability condition)

frame_diff    = diff(analysis_frames);
frame_beg     = find(frame_diff  >  0.1) +1;
frame_end     = find(frame_diff  < -0.1);
% redefine the begin and end index of frames
N_frames      = length(frame_end);
%#############################################################
%############################################################# cutting
% In case where analysis frame exceeds max_frame is cut

for i=1:N_frames
   
    t1 = frame_beg(i);
    t2 = frame_end(i);
        
    while t2-t1+1  > max_frame        
     
        t1 = t1 +  max_frame + n;        
        analysis_frames(t1-n:t1-1) = 0;        
        d_f(t1-n:t1-1)     = 0;
        d_b(t1-n:t1-1)     = 0;
        d_f_org(t1-n:t1-1) = 0;
        d_b_org(t1-n:t1-1) = 0;
        display('max frame detected')
    end
    
end

frame_diff    = diff(analysis_frames);
frame_beg     = find(frame_diff  >  0.1) +1;
frame_end     = find(frame_diff  < -0.1);

N_frames      = length(frame_end);

%#############################################################
%############################################################# combining

pattern_counter       = zeros(10,1); 

N = length(d_f);

d_local_fusion_rules        = zeros(N,1);
d_global_union_rules        = zeros(N,1);
d_global_intersecion_rules  = zeros(N,1);


for i=1:N_frames
    
    
    %###########################################
    frame_f = d_f(frame_beg(i)-1:frame_end(i)+1);
    frame_b = d_b(frame_beg(i)-1:frame_end(i)+1);
    
    frame_f_diff = diff(frame_f);
    alarm_f_beg  = find(frame_f_diff  >  0.1 ) +1;
    alarm_f_end  = find(frame_f_diff  < -0.1 );
    
    
    frame_b_diff = diff(frame_b);
    alarm_b_beg  = find(frame_b_diff  >  0.1 ) +1;
    alarm_b_end  = find(frame_b_diff  < -0.1 );
        
    
    alarms_f = length(alarm_f_end); %the numbers of forwards alarm
    alarms_b = length(alarm_b_end);% the numbers of backwards alarm
    %###########################################
    
   
    %###########################################    
    if isempty(alarm_f_beg) == false
        t1_f = frame_beg(i) + alarm_f_beg(1)    - 2;
        t2_f = frame_beg(i) + alarm_f_end(end)  - 2;
    else
        t1_f = 0;
        t2_f = 0;
    end

    if isempty(alarm_b_beg) == false
        t1_b = frame_beg(i) + alarm_b_beg(1)    - 2;
        t2_b = frame_beg(i) + alarm_b_end(end)  - 2;                    
    else
        t1_b = 0;
        t2_b = 0;
    end
    %###########################################           
                
    if alarms_f <= 1 && alarms_b <= 1 %elementary detection patterns
        
                      
                %############################# subclass A1
                if t1_b == t1_f && t2_b == t2_f
                   
                    pattern = 1;
                   
                    t1_A = frame_beg(i); 
                    t2_A = frame_end(i);
                                        
                end
                
                %############################# subclass A2
                if (t1_b >= t1_f && t2_b < t2_f) || (t1_b > t1_f && t2_b <= t2_f)
                   
                    pattern = 2; 
                          
                    t1_A = t1_f; 
                    t2_A = t2_b;
                end
                
                %############################# subclass A3
                if (t1_b <= t1_f && t2_b > t2_f) || (t1_b < t1_f && t2_b >= t2_f)
                   
                    pattern = 3; 
                    
                    t1_A = t1_f; 
                    t2_A = t2_b; 
                end
                
                %############################# subclass A4
                if t1_b > t1_f && t1_b <= t2_f && t2_b > t2_f
                   
                    pattern = 4; 
                    
                    t1_A = t1_f; 
                    t2_A = t2_b; 
                end
                
                %############################# subclass A5
                if t1_f > t1_b && t1_f <= t2_b && t2_f > t2_b
                   
                    pattern = 5; 
                    
                    t1_A = t1_f; 
                    t2_A = t2_b;  
                end
                
                %############################# subclass B1
                if t1_b > t2_f && t1_b - t2_f <= n
                   
                    pattern = 6; 
                   
                   t1_A = t1_f; 
                   t2_A = t2_b; 
                end
                
                %############################# subclass B2
                if t1_f > t2_b && t1_f - t2_b <= n
                  
                    pattern = 7; 
                    
                    t1_A = t1_b; 
                    t2_A = t2_f; 
                end
                
                %############################# subclass C1
                if t1_f > 0 && t2_f > 0  && t1_b == 0 && t2_b == 0
                   
                    pattern = 8; 
                                       
                    %######################################################
                    frame_f_org         = d_f_org(frame_beg(i)-1:frame_end(i)+1);                    
                    frame_int           = frame_f & frame_f_org;
                    frame_int_diff      = diff(frame_int);
                    alarm_f_beg_org     = find(frame_int_diff  >  0.1) + 1;  
                    
                    t1_org = frame_beg(i) + alarm_f_beg_org - 2;
                    
                    extension = 0; 
                    for j=1:delta        
                        if i+1 <= N_frames  
                            if frame_beg(i+1)- t1_org - j > n
                               extension = j;
                            end
                        else
                            extension = delta;
                        end
                    end
                    %######################################################
                    
                    
                    t1_A = t1_f; 
                    t2_A = t1_org + extension; 
                    
                end
                
                %############################# subclass C2
                if t1_f == 0 && t2_f == 0  && t1_b > 0 && t2_b > 0
                  
                    pattern = 9; 
                    
                    %######################################################
                    frame_b_org         = d_b_org(frame_beg(i)-1:frame_end(i)+1);                    
                    frame_int           = frame_b & frame_b_org;
                    frame_int_diff      = diff(frame_int);
                    alarm_b_end_org     = find(frame_int_diff  <  -0.1);  
                    
                    t2_org = frame_beg(i) + alarm_b_end_org - 2;
                    
                    extension = 0; 
                    for j=1:delta        
                        if i-1 > 0  
                            if  t2_org - j - frame_end(i-1) > n
                               extension = j;
                            end
                        else
                            extension = delta;
                        end
                    end
                    
                    t1_A = t2_org - extension; 
                    t2_A = t2_b;
                    
                end
                %#############################
    
    else
                %#############################
                %complex patterns
                pattern = 10; 
                
                t1_A = t1_f; 
                t2_A = t2_b ; 
                %#############################
    end
    

    d_local_fusion_rules(t1_A:t2_A)                           = 1;
    d_global_union_rules(frame_beg(i):frame_end(i))           = 1;    
    d_global_intersecion_rules(frame_beg(i)-1:frame_end(i)+1) = (frame_f & frame_b);            
    
    pattern_counter(pattern) = pattern_counter(pattern) + 1;
 
    
%###########################################################  
if 1<0
% The figure shows the result of combining forward-time and backward-time detection alarms

offset = 2*n;

hold off
stairs(d_f_org(frame_beg(i)-offset:frame_end(i)+offset)+10)
hold on
stairs(d_b_org(frame_beg(i)-offset:frame_end(i)+offset)+8)
hold on
stairs(d_f(frame_beg(i)-offset:frame_end(i)+offset)+6)
hold on
stairs(d_b(frame_beg(i)-offset:frame_end(i)+offset)+4)
hold on
stairs(d_local_fusion_rules(frame_beg(i)-offset:frame_end(i)+offset))

legend('d_f','d_b','d_f^*','d_b^*','d_{fb}')
pause

end
%###########################################################   
    
end

%##########################################################################




end