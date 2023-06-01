clear all;close all;clc;
in_directory= './music_dataset';
out_directory = './music_dataset';

[impul_noise, Fs] =  audioread([in_directory '/NoiseFromPaper.wav']);

for rec = 1: 20
    filename = [in_directory '/rec_ORG_' num2str(rec) '.wav'];
    org = audioread(filename);
    ref = org + impul_noise;
    outname = [out_directory '/rec_REF_' num2str(rec) '.wav'];
    audiowrite(outname, ref, Fs);
end

