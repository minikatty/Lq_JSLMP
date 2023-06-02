clear;close all;close all;
% addpath('D:\Program Files\MATLAB\R2016b\JSTSP\speech_dataset');
%% Algorithm
clc;
ov=4;                                                      %overlap factor
inc=128;                                                 %increment
nw=inc*ov;                                             %window length
W=hamming(nw,'periodic');                %hamming window
W=W/sqrt(sum(W(1:inc:nw).^2));        %normalize window                 
method='spJSLMP';%'spJSLMP';'LpADM';
% SNR=[0 5 10 15];
% dir = ['.\speech_out\female\' num2str(SNR(4)) 'db\corrupted\'];
% tic
in_indirectory = './dataset/';
for num =1:20    
    filename = [in_indirectory 'bgn_REF_' num2str(num) '.wav']; %'rec_REF_'
%     filename2 = [in_indirectory 'rec_ORG_' num2str(num) '.wav'];
%     filename3 = [in_indirectory 'rec_REF_' num2str(num) '.wav'];
%     [ORG,~]=audioread(filename2);
    [SAR,fs]=audioread(filename);
%     speechfile=['speech_dataset/male/10db/','Ms_REF',num2str(num),'_10dB.wav'];
%     [s,fs]=audioread(speechfile);
%     [REF,~]=audioread(filename3);
%     filename2=['rec_ORG_',num2str(i),'.wav'];
%     [clean,~]=audioread(filename2); 
%     ORG_frame = enframe(ORG,W,inc);
%     REF_frame = enframe(REF,W,inc);
%     SAR_frame=enframe(SAR,W,inc);
    Y=enframe(SAR,W,inc);%REF
%    Y=enframe(s,W,inc); %enframe,specify window type and inc,ref
    rec=zeros(size(Y)); 
%     rec2=rec;
    p=1.6;q1=1;q2=q1;
%      tic
    for j=1:size(Y,1)
%         [count,~]=DetectImpul(Y(j,:));
        if strcmp(method,'spJSLMP')
            if j==1
                [rec(j,:),flag]=spJSLMP(Y(j,:),q1,q2,p,'normal');                 %JSLMP DCT filter 
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
%       toc
    
X=v_overlapadd(rec,W,inc);              %reconstruct
%     filename3=['rec_raw_',num2str(i),'.wav'];
    filename3=['.\JSLMP4Data2\bgn_rec_jslmp_',num2str(num),'.wav'];%jslmp171
    audiowrite(filename3,X,fs);
end
% synthesis;