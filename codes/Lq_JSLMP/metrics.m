
%% Add noise
%%  Calculate the STOI of noise audio
clc;clear;close;
num=10;
SNR=[0 5 10 15];
path=['.\speech_out\female\' num2str(SNR(4)) 'db\'];
for i=1:10
    orgfile=[path,'org\Fs_org_red',num2str(i),'.wav'];
    [org,fs]=audioread(orgfile);
    recfile=[path,'jslmp\Fs_rec_jslmp_' num2str(i) '_' num2str(SNR(4)) 'dB.wav'];
    %[path,'jslmp171\Ms_rec_',num2str(i),'_15dB.wav'];
    %[path,'corrupted\Ms_REF',num2str(i),'_15dB.wav'];
    %[path,'lpADM\Ms_rec_lpAdm_',num2str(i),'_15dB.wav'];
    %
    [est,~]=audioread(recfile);
    sSNR(i)=v_snrseg(est,org,fs);
    STOI(i)=stoi(org, est, fs);
    SDR(i)=eval_sdr2(est,org,fs);
end
mean(STOI)