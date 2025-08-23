clear;close all;close all;
% addpath('D:\Program Files\MATLAB\R2016b\JSTSP\speech_dataset');
%% Algorithm
clc;
ov=2;                                                      %overlap factor 4
inc=256;                                                 %increment
nw=inc*ov;                                             %window length
W=hamming(nw,'periodic');                %hamming window
W=W/sqrt(sum(W(1:inc:nw).^2));        %normalize window                 
method='spJSLMP';%'spJSLMP';'LpADM';
% SNR=[0 5 10 15];
% dir = ['.\speech_out\female\' num2str(SNR(4)) 'db\corrupted\'];
% tic

% in_indirectory = 'noisy_dataset';
% out_indirectory = 'Lq_JSLMP_speech_results_val';
out_indirectory = [];
% clean_dir = fullfile(pwd, 'clean_speech');
% if ~exist(out_indirectory, 'dir')
%     mkdir(out_indirectory);
% end
% files = dir(fullfile(in_indirectory, '*_0dB.wav'));
fs = 16e3;
% frame_len = 512;
% frame_num = 61;% 1s-duration data
prob = 0.2;
y1 = 0.01*randn(fs,1); % test memery
y2 = 0.15*randn(fs,1);

num_y2 = round(prob * fs);   % 来自 y2 的帧数
num_y1 = fs - num_y2;       % 来自 y1 的帧数

% 从 y1 和 y2 中随机抽样对应数量的帧
idx_y1 = randperm(fs, num_y1);
idx_y2 = randperm(fs, num_y2);

y_mix = [y1(idx_y1); y2(idx_y2)];
y_mix = y_mix(randperm(length(y_mix)));  % 打乱混合顺序
%##########################################################################
%##########################################################################
% y = y_mix(:);
for num =1:1
    %length(files)
    % filename = files(num).name; %'rec_REF_'
%     filename2 = [in_indirectory 'rec_ORG_' num2str(num) '.wav'];
%     filename3 = [in_indirectory 'rec_REF_' num2str(num) '.wav'];
%     [ORG,~]=audioread(filename2);
    % path = fullfile(in_indirectory, filename);
    data = y_mix;
    % [data,fs]=audioread(path); % column vector
%     speechfile=['speech_dataset/male/10db/','Ms_REF',num2str(num),'_10dB.wav'];
%     [s,fs]=audioread(speechfile);
%     [REF,~]=audioread(filename3);
%     filename2=['rec_ORG_',num2str(i),'.wav'];
%     [clean,~]=audioread(filename2); 
%     ORG_frame = enframe(ORG,W,inc);
%     REF_frame = enframe(REF,W,inc);
%     SAR_frame=enframe(SAR,W,inc);
    Y=enframe(data,W,inc);%REF
%    Y=enframe(s,W,inc); %enframe,specify window type and inc,ref
    rec=zeros(size(Y)); 
%     rec2=rec;
    p=1;q1=0.7;q2=q1;
%      tic
    tic
    for j=1:size(Y,1)
%         [count,~]=DetectImpul(Y(j,:));
        if strcmp(method,'spJSLMP')
            if j==1
                [rec(j,:),flag]=spJSLMP(Y(j,:),q1,q2,p,'normal');
                %JSLMP DCT filter 
%                 [rec2(j,:),flag2]=Pro_spJSLMP(Y(j,:),q1,q2,p,'normal'); 
            else
                [rec(j,:),flag]=spJSLMP(Y(j,:),q1,q2,p,'normal',rec(j-1,:)');%
%                 [rec2(j,:),flag2]=Pro_spJSLMP(Y(j,:),q1,q2,p,'normal'); 
            end
        elseif strcmp(method,'LpADM')
            if j==1
                [rec(j,:)]=LpADM(Y(j,:)',0.9,1);                 
            else
                rec(j,:)=LpADM(Y(j,:)',0.9,1,rec(j-1,:)');
            end
            
        end
          
%           title('Origina');
%           proposed = rec(j,:);  SAR_re = SAR_frame(j,:);
%           org = ORG_frame(j,:); ref = REF_frame(j,:);
%           figure;plot(org,'r');xlim([0 length(rec(j,:))]);hold on;
%           plot(proposed);xlim([0 length(proposed)]);hold on;
%           plot(SAR_re);xlim([0 length(rec(j,:))]);
%           hold on;
%           plot(ref,'-');xlim([0 length(rec(j,:))]);
%          hold on;plot(rec(j,:));xlim([0 length(rec(j,:))]);
%          legend('ORG','Pro','SAR','Raw');
%          error1=sum((rec(j,:)-Y1(j,:)).^2);% proposed error
%          title('Proposed_rec');
%           figure;plot(Y2(j,:),'k');xlim([0 length(rec(j,:))]);title('GT');
%           hold  on;plot(Y(j,:),'r');
         
%          error0=sum((ref-org).^2);% raw errror
%          error1=sum((proposed-org).^2);% jsLMP
%          error2=sum((SAR_re-org).^2);%SAR
         
%          [error1 error2] %REF proposed SAR           
%           figure;
%           plot(Y(j,:));xlim([0 length(rec(j,:))]);hold on;
%           plot(Y1(j,:),'r');legend('SAR','JS','Ref','true');
    end
    runtime = toc;
    runtime/size(Y,1)*1000
%       toc    
    X=v_overlapadd(rec,W,inc);              %reconstruct
    deg = X;
    %     filename3=['rec_raw_',num2str(i),'.wav'];
    % output_name = strrep(filename, 'noisy', 'rec');
    output_name = 'test_perf.wav';
    out_path = fullfile(out_indirectory, output_name);
    % filename3=['.\JSLMP4Data2\bgn_rec_jslmp_',num2str(num),'.wav'];%jslmp171
    audiowrite(out_path,X,fs);
    %     % 判定干净语音的路径
    % if contains(filename, '_red_rec_') && ismember(snr, [0, 5])
    %     %处理后的语音是_red_rec_
    %     %原始带噪语音室_red_noisy_
    %     clean_name = sprintf('%ss_org_red_again_%s.wav', filename(1), utt_index);
    % else
    %     clean_name = sprintf('%ss_org_red%s.wav', filename(1), utt_index);
    % end
    % clean_path = fullfile(clean_dir, clean_name);
    % ref = audioread(clean_path,X,fs);
    % PESQ = py.pesq(fs, ref, deg, 'wb');
    
end
% synthesis;