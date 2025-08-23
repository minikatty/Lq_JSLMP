import os
import csv
from pesq import pesq
from scipy.io import wavfile
from collections import defaultdict
import re

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
clean_dir = os.path.join(BASE_DIR,"clean_speech")
BASE_DIR = os.path.dirname(BASE_DIR)
algorithm = 'Lq_JSLMP2' #SAR, Lq_JSLMP, LpADM，noisy,gtcrn, ft-gtcrn
noise_source = 'measured'# measured,[],paper
enhanced_dir = os.path.join(BASE_DIR,"Lq_JSLMP_speech_results_val")
# 'noisy_dataset','Lq_JSLMP_speech_results','eh_noisy_dataset_results_withFinetuning'
# 'eh_noisy_dataset_results_withoutFinetuning','LpADM_speech_results','SAR_speech_results'
#eh_noisy_dataset_results_withoutFinetuning，
# eh_noisy_dataset_results_withFinetuning
output_csv = 'val_'+algorithm +'_' + noise_source +'_PESQ.csv'

# output_csv = "eh_noFineTuning_pesq_scores_paper_red05.csv"# eh_noFineTuning_pesq_scores_paper_red05.csv

results = []
snr_groups = defaultdict(list)  # 用于按 SNR 分组
# for no scale data
"""
for fname in os.listdir(enhanced_dir):
    if not fname.endswith(".wav") :
        continue
    if "_rec" not in fname:
        continue
    gender = fname[0]  # F or M
    index_part = fname.split("_")[1]  # e.g., rec1
    utt_index = index_part[6:]  # e.g., '1'
"""
# for scaled data:0,5dB paper noise

for fname in os.listdir(enhanced_dir):
    if not (f"{noise_source}_rec_" in fname or f"{noise_source}_red_rec_" in fname):
        continue
    # _noisy_,_rec_
    # if "_red" not in fname:
    #     continue
    gender = fname[0]  # F or M 有点多余
    index_part = fname.split("_")[1]  # e.g., rec1
    utt_index = index_part[6:]  # e.g., '1'
    """
        gender = fname[0]  # F or M
        index_part = fname.split("_")[1]  # e.g., rec1
        utt_index = index_part[3:]  # e.g., '1'
    """
    # 提取 SNR
    match = re.search(r'rec_(\d+)dB', fname)
    if not match:
        print(f"无法提取 SNR：{fname}")
        continue
    snr = int(match.group(1))

    # 构造参考文件名，如 Fs_org_red1.wav
    if "_red_rec" in fname and snr in [0,5]:
        clean_name = f"{gender}s_org_red_again_{utt_index}.wav"
    else:
        clean_name = f"{gender}s_org_red{utt_index}.wav"

    clean_path = os.path.join(clean_dir, clean_name)
    enh_path = os.path.join(enhanced_dir, fname)

    if not os.path.exists(clean_path):
        print(f"缺失干净语音：{clean_path}")
        continue

    fs_ref, ref = wavfile.read(clean_path)
    fs_deg, deg = wavfile.read(enh_path)

    if fs_ref != 16000 or fs_deg != 16000:
        print(f"非16kHz采样率：{fname}")
        continue

    try:
        score = pesq(fs_ref, ref, deg, 'wb')
        results.append([fname, score])  # 保存 snr 信息
        snr_groups[snr].append(score)  # 分组累积
        print(f"{fname} → PESQ = {score:.4f}")
    except Exception as e:
        print(f"计算失败：{fname}, 原因：{e}")

# 逐条写入 CSV 文件
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["filename", "pesq"])
    writer.writerows(results)

# 写入对应SNR音频的pesq均值
avg_output_csv = 'eh_'+algorithm+'_'+noise_source+'_avg_PESQ.csv'
# "eh_noFineTuning_pesq_avg_per_snr.csv"
with open(avg_output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["snr(dB)", "avg_pesq"])
    for snr in sorted(snr_groups.keys()):
        avg_score = sum(snr_groups[snr]) / len(snr_groups[snr])
        writer.writerow([snr, round(avg_score, 4)])
        print(f"SNR {snr} dB → 平均 PESQ = {avg_score:.4f}")

print(f"\n所有单个文件的 PESQ 分数已保存到 {output_csv}")
print(f"各 SNR 平均 PESQ 分数已保存到 {avg_output_csv}")

print(f"\n所有 PESQ 分数已保存到 {output_csv}")
