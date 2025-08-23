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
% algorithm = 4 - the semi-causal, equipped with decision-feedback detection scheme,
% and 4 bidirectional algorithms i.e., theirs noncausal extensions.
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
% y_f  - the processed signal in forward-time
% y_fb - the processed signal, bidirectional signal processing
%
% Date: 08.2016
% ? All rights reserved.
%##########################################################################
%##########################################################################
clear all;close all;clc
% input_directory     = ['./noisy_dataset/'];
output_directory    = ['./SAR_speech_results/'];
% 
% files = dir(fullfile(input_directory, '*.wav'));
% 
% for rec =1:length(files)
%     fname = files(rec).name;
%     if ~contains(fname,'_15dB')
%     % _paper_:论文中提取出俩的
%     % _measured_:我们实际录到的
%     %确认处理音频是恢复出来的文件,因为还有noisy_data要处理
%     % 带噪数据用_noisy_
%         continue;
%     end
% rec         = 10;
% file        = ['Archive audio ' num2str(rec)];
% file   = files.name;%['bgn_REF_' num2str(rec)];
% path = fullfile(input_directory, fname);
%[input_directory file '.wav'];
% [track,Fs]  = audioread(path);
% fname = 'test_perfor.wav';
Fs= 16000;
% y           = track(:,1);
frame_len = 512;
frame_num = 61;% 1s-duration data
p = 0.2;
y1 = 0.01*randn(frame_num,frame_len); % test memery
y2 = 0.15*randn(frame_num,frame_len);

num_y2 = round(p * frame_len);   % 来自 y2 的帧数
num_y1 = frame_len - num_y2;       % 来自 y1 的帧数

% 从 y1 和 y2 中随机抽样对应数量的帧
idx_y1 = randperm(frame_len, num_y1);
idx_y2 = randperm(frame_len, num_y2);

y_mix = [y1(:, idx_y1), y2(:, idx_y2)];
y_mix = y_mix(:, randperm(size(y_mix,2)));  % 打乱混合顺序
%##########################################################################
%##########################################################################
y = y_mix(:);

%##########################################################################
%##########################################################################

% Simulation settings:
k         = 240;
M         = 2*k+1; % segment length
n         = 10;    % the order of  AR model


mu_alpha  = 5;   % prediction based detection threshold
mu_beta   = 4.5;   % interpolation based detection threshold mu_beta >=mu_alpha
kmax      = 60;    % maximum length of detection alarm

threshold_alpha = mu_alpha^2;
threshold_beta  = mu_beta^2;

algorithm = 4;

%##########################################################################
%##########################################################################
tic            
[y_f,y_fb,d_fb] = bidirectional_signal_processing(y,n,M,threshold_alpha,threshold_beta,kmax,algorithm);
% index_set_BGN(rec,:) = d_fb'; % placed by row
% Ftest1 = [output_directory  file '_' num2str(algorithm) '_' num2str(mu_alpha)   '_f.wav'];
% audiowrite(Ftest1,y_f,Fs); 
% F_REF   = [input_directory  'bgn_REF_' num2str(rec) '.wav'];
% file1 = ['bgn_rec_' num2str(rec) '_'  num2str(algorithm) '_' num2str(mu_alpha)   '_fb.wav'];
% output_name = strrep(fname, 'noisy', 'rec');
output_name = 'test_perfor.wav';
Ftest2 = [output_directory  output_name];%'_' num2str(algorithm) '_' num2str(mu_alpha)   '_fb.wav'
audiowrite(Ftest2,y_fb,Fs);
runtime = toc;
runtime/frame_num*1000
% end    
if 1<0 %PEAQ tool, only when the orignal (clean) recording is available


F_ORG   = [input_directory  'rec_ORG_' num2str(rec) '.wav'];

StartS = Fs+1;
EndS   = 21*Fs;

ODG1   = PQevalAudio (F_ORG, F_REF, StartS, EndS);
ODG2   = PQevalAudio (F_ORG, Ftest2, StartS, EndS);

ODG    = [ODG1 ODG2]

save([output_directory 'ODG_restored_' num2str(rec) '_' num2str(algorithm) '_' num2str(mu_alpha) '.mat'],'ODG')

end

%##########################################################################
%##########################################################################
if 1<0 %PEAQ tool, only when orignal (clean) recording is available

    display('Ground truth reconstruction')
    load([input_directory 'noise_pulse_localization.mat']) % d_org
    
    F_REF   = [input_directory  'rec_REF_' num2str(rec) '.wav'];
    [y,Fs]  = audioread(F_REF); 
    [y_gt]  = reconstruction_algorithm(y,n,M,d_org);
    
    Ftest2 = [input_directory  'rec_GT_' num2str(rec) '.wav'];
    audiowrite(Ftest2,y_gt,Fs); 
    
    F_ORG   = [input_directory  'rec_ORG_' num2str(rec) '.wav'];
    F_REF   = [input_directory  'rec_REF_' num2str(rec) '.wav'];
    F_GT    = [input_directory  'rec_GT_'  num2str(rec) '.wav'];    
    
    
    StartS  = Fs+1;% StartS  = Fs+1;
    EndS    = 21*Fs;% EndS    = 21*Fs;
    
    ODG1    = PQevalAudio (F_ORG, F_GT , StartS, EndS);
    ODG2    = PQevalAudio (F_ORG, F_REF, StartS, EndS);
    
    ODG     = [ODG2 ODG1 ]; %    
    
    save([output_directory 'ODG_ground_truth_' num2str(rec) '.mat'],'ODG')
end
%##########################################################################
%##########################################################################
% end
