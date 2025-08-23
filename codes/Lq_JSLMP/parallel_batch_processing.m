% 多核并行版本 - 批量运行 spJSLMP 或 LpADM
clear; clc; close all;
delete(gcp('nocreate'));
% 设置并行核数
nCores = 8; % 可修改核心数
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local', nCores);
end
%% 参数设置
ov = 4;                      % overlap factor
inc = 128;                   % increment
nw = inc * ov;              % window length
W = hamming(nw, 'periodic'); 
W = W / sqrt(sum(W(1:inc:nw).^2));  % normalize window
method = 'spJSLMP';          % 可选: 'spJSLMP' 或 'LpADM'
p = 1.6; q1 = 0.7; q2 = q1;



%% 路径设置
in_indirectory = 'noisy_dataset';
out_indirectory = ['Lq_JSLMP_speech_' ...
    'results_V2'];
if ~exist(out_indirectory, 'dir')
    mkdir(out_indirectory);
end
files = dir(fullfile(in_indirectory, '*_15dB.wav'));

%% 并行处理
parfor num = 1:length(files)
    filename = files(num).name;
    input_path = fullfile(in_indirectory, filename);
    [data, fs] = audioread(input_path);

    % 分帧
    Y = enframe(data, W, inc);
    rec = zeros(size(Y));

    % 帧级处理
    for j = 1:size(Y,1)
        if strcmp(method,'spJSLMP')
            if j == 1
                [rec(j,:), ~] = spJSLMP(Y(j,:), q1, q2, p, 'normal');
            else
                [rec(j,:), ~] = spJSLMP(Y(j,:), q1, q2, p, 'normal', rec(j-1,:)');
            end
        else strcmp(method,'LpADM')
            if j == 1
                rec(j,:) = LpADM(Y(j,:)', 0.9, 1);
            else
                rec(j,:) = LpADM(Y(j,:)', 0.9, 1, rec(j-1,:)');
            end
        end
    end

    % 重构并保存
    X = v_overlapadd(rec, W, inc);
    output_name = strrep(filename, 'noisy', 'rec');
    out_path = fullfile(out_indirectory, output_name);
    audiowrite(out_path, X, fs);
    fprintf('完成: %s\n', filename);
end

disp('全部处理完成');