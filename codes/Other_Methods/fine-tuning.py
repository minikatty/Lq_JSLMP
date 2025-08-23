import os
import random
import torch
import torch.nn as nn
import torch.optim as optim
import soundfile as sf
from tqdm import tqdm
from gtcrn import GTCRN
from loss import HybridLoss
from torch.utils.data import DataLoader,Dataset
from torch.optim.lr_scheduler import ReduceLROnPlateau

test_results_dir = 'testdata_results_with_fine_tuning' #results_without_fine_tuning
os.makedirs(test_results_dir, exist_ok=True)

device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
"""
输入类型就是按（B, T, F, 实虚部）
"""
class SpeechDataset(Dataset):
    def __init__(self, noisy_dir, clean_dir, file_list, sample_rate, n_fft, hop_length, win_length, window):
        self.noisy_dir = noisy_dir
        self.clean_dir = clean_dir
        self.file_list = file_list
        self.sample_rate = sample_rate
        self.n_fft = n_fft
        self.hop_length = hop_length
        self.win_length = win_length
        self.window = window

    def __len__(self):
        return len(self.file_list)

    def __getitem__(self, idx):
        fname = self.file_list[idx]
        # e.g. Fs_REF10_5dB.wav -> Fs_org_red10.wav
        name_parts = os.path.splitext(fname)[0].split("_")
        speaker = name_parts[0]
        index = name_parts[1][3:]
        clean_name = f"{speaker}_org_red{index}.wav"

        noisy_path = os.path.join(self.noisy_dir, fname)
        clean_path = os.path.join(self.clean_dir, clean_name)

        waveform_noisy, fs = sf.read(noisy_path)
        waveform_clean, fs2 = sf.read(clean_path)
        assert fs == self.sample_rate and fs2 == self.sample_rate

        stft_noisy = torch.stft(torch.from_numpy(waveform_noisy).to(device).float(),
                                self.n_fft, self.hop_length, self.win_length,
                                window=self.window, return_complex=False)

        stft_clean = torch.stft(torch.from_numpy(waveform_clean).to(device).float(),
                                self.n_fft, self.hop_length, self.win_length,
                                window=self.window, return_complex=False)

        return stft_noisy, stft_clean


# ---------- 参数设置 ----------
noisy_dir = 'test_data_without_fine_tuning'
clean_dir = 'clean_speech'
checkpoint_path = '../checkpoints/model_trained_on_dns3.tar'
output_dir = '../checkpoints'
save_name = 'model_finetuned_on_dns3_v2.pth'
epochs = 5
batch_size = 2  # 由于变长语音，设为1
lr = 1e-4
sample_rate = 16000

# ---------- STFT ----------
n_fft = 512
hop_length = 256
win_length = 512
window = torch.hann_window(win_length).pow(0.5).to(device)

# ---------- train/test dataset ----------
all_files = sorted([f for f in os.listdir(noisy_dir) if f.endswith('.wav')])
random.seed(42)
random.shuffle(all_files)
train_files = all_files[:60] # 60 for fine-tuning, 20 for validation
test_files = all_files[60:]

# ------------ save test dataset for evaluation------------
test_list_path = os.path.join(test_results_dir, 'test_files_used_for_finetuning_eval.txt')
with open(test_list_path, 'w') as f:
    for fname in test_files:
        f.write(fname + '\n')
print(f"测试文件列表已保存至: {test_list_path}")

train_dataset = SpeechDataset(noisy_dir, clean_dir, train_files, sample_rate, n_fft, hop_length, win_length, window)
test_dataset  = SpeechDataset(noisy_dir, clean_dir, test_files,  sample_rate, n_fft, hop_length, win_length, window)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
test_loader  = DataLoader(test_dataset,  batch_size=batch_size, shuffle=False)

# ---------- load trained model ----------

model = GTCRN().to(device)
ckpt = torch.load(checkpoint_path, map_location=device)
model.load_state_dict(ckpt['model'])
model.train()

# ------- Hyperbird Loss & optimizer -------------
loss_fn = HybridLoss().to(device)
optimizer = optim.Adam(model.parameters(), lr=lr)
scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=3, verbose=True)

# ---------- 训练循环 ----------
for epoch in range(epochs):
    total_loss = 0.0
    model.train()
    for noisy_stft, clean_stft in tqdm(train_loader, desc=f"Epoch {epoch+1}/{epochs}"):
        noisy_stft = noisy_stft.to(device)
        clean_stft = clean_stft.to(device)
        pred_stft = model(noisy_stft) # 输出维度: (B, F, T, 2)

        loss = loss_fn(pred_stft, clean_stft)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    avg_loss = total_loss / len(train_files)
    print(f"Epoch {epoch+1} average loss: {avg_loss:.6f}")

    # ---------- validation ----------
    model.eval()
    with torch.no_grad():
        val_loss = 0.0
        for fname in test_files:
            for noisy_stft, clean_stft in test_loader:
                noisy_stft = noisy_stft.to(device)
                clean_stft = clean_stft.to(device)

                pred_stft = model(noisy_stft)

                loss = loss_fn(pred_stft, clean_stft)
                val_loss += loss.item()
        val_loss /= len(test_loader)
        scheduler.step(val_loss)
        print(f"Validation loss: {val_loss:.6f}")

# ---------- save_model ----------
torch.save(model.state_dict(), os.path.join(output_dir, save_name))
print(f"模型已保存至 {save_name}")
