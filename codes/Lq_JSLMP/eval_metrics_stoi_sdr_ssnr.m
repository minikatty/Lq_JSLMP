clc;clear;close;
% 评估 STOI、SDR、SSNR，并按 SNR 分组平均
clean_dir = fullfile(pwd, 'clean_speech');
noise_source = 'measured';% measured,[],paper

enhanced_dir = fullfile(pwd, 'eh_noisy_dataset_results_withFinetuning');
% pwd:当前文件夹,'Lq_JSLMP_speech_results_v2'
% ‘noisy_dataset’,'Lq_JSLMP_speech_results','eh_noisy_dataset_results_withFinetuning'
% 'eh_noisy_dataset_results_withoutFinetuning','LpADM_speech_results','SAR_speech_results'
algorithm = 'ft-gtcrn';%SAR, Lq_JSLMP, LpADM，noisy,gtcrn, ft-gtcrn,Lq_JSLMP_v2
output_csv = ['eh_' algorithm '_' noise_source '_metrics.csv'];%
avg_csv = ['eh_' algorithm '_' noise_source '_avg_metrics.csv'];

files = dir(fullfile(enhanced_dir, '*.wav'));
results = {};
snr_map = containers.Map('KeyType','double','ValueType','any');

for i = 1:length(files)
    fname = files(i).name;
    if ~contains(fname, [noise_source '_rec_'])&&~contains(fname, [noise_source '_red_rec_'])
    % _paper_:论文中提取出的
    % _measured_:我们实际录到的
    %确认处理音频是恢复出来的文件,因为还有noisy_data要处理
    % 带噪数据用_noisy_
        continue;
    end

    % 提取 SNR 值
    snr_expr = regexp(fname, '_rec_(\d+)dB', 'tokens');%带噪语音为_noisy_,语音增强的为_rec_
    if isempty(snr_expr)
        fprintf('无法提取 SNR：%s\n', fname);
        continue;
    end
    snr = str2double(snr_expr{1}{1});

    % 提取语音编号
    utt_index_expr = regexp(fname, '_speech(\d+)_', 'tokens');
    utt_index = utt_index_expr{1}{1};
    % gender = fname(1);  % 提取说话人性别F or M

    % 判定干净语音的路径
    if contains(fname, '_red_rec_') && ismember(snr, [0, 5])
        %处理后的语音是_red_rec_
        %原始带噪语音室_red_noisy_
        clean_name = sprintf('%ss_org_red_again_%s.wav', fname(1), utt_index);
    else
        clean_name = sprintf('%ss_org_red%s.wav', fname(1), utt_index);
    end

    clean_path = fullfile(clean_dir, clean_name);
    enh_path = fullfile(enhanced_dir, fname);

    if ~isfile(clean_path)
        fprintf('缺失干净语音：%s\n', clean_path);
        continue;
    end

    [ref, fs1] = audioread(clean_path);
    [deg, fs2] = audioread(enh_path);

    if fs1 ~= 16000 || fs2 ~= 16000
        fprintf('非16kHz采样率：%s\n', fname);
        continue;
    end

    % 计算指标
    try
        s_stoi = stoi(ref, deg, fs1); % STOI(i)=stoi(org, est, fs);
        s_sdr = eval_sdr2(deg,ref,fs1); % SDR(i)=eval_sdr2(est,org,fs);
        s_ssnr = v_snrseg(deg,ref,fs1); % v_snrseg(est,org,fs);
        s_pesq = py.pesq.pesq(int32(fs1), py.numpy.array(ref), py.numpy.array(deg), 'wb');
        results(end+1, :) = {fname, snr, s_pesq, s_stoi, s_sdr, s_ssnr};

        if ~isKey(snr_map, snr)
            snr_map(snr) = [];
        end

        snr_map(snr) = [snr_map(snr); s_pesq, s_stoi, s_sdr, s_ssnr];
        fprintf(['Clean speech: %s, Enhanced speech:%s\n'...
            'PESQ =%.4f, STOI=%.4f, SDR=%.4f, SSNR=%.4f\n'], ...
            clean_name, fname, s_pesq, s_stoi, s_sdr, s_ssnr);
    catch ME
        fprintf('计算失败：%s，原因：%s\n', fname, ME.message);
    end
end

% 写入每一条结果
fid = fopen(output_csv, 'w');
fprintf(fid, 'filename,snr,pesq,stoi,sdr,ssnr\n');
for i = 1:size(results,1)
    fprintf(fid, '%s, %d, %.4f, %.4f, %.2f, %.2f\n', results{i,1}, ...
        results{i,2}, results{i,3}, results{i,4}, results{i,5}, results{i,6});
end
fclose(fid);

% 写入均值
fid = fopen(avg_csv, 'w');
fprintf(fid, 'snr,avg_pesq,avg_stoi,avg_sdr,avg_ssnr\n');
keys = snr_map.keys;
for i = 1:length(keys)
    snr = keys{i};
    data = snr_map(snr);
    avg_vals = mean(data,1);
    fprintf(fid, '%d,%.4f,%.4f,%.2f,%.2f\n', snr, ...
        avg_vals(1), avg_vals(2), avg_vals(3), avg_vals(4));
end
fclose(fid);