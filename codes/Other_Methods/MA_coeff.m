clc;
clear;
close all;

% 生成一个简单的噪声信号
n = 1000; % 信号长度
% v = sqrt(0.1)*randn(n,1); % 生成白噪声
order_ma = 4; % 假设MA模型阶数为4
%% 产生均值为0方差为0.1的噪声信号
v = sqrt(0.003)*randn(n,1);

%% 产生AR模型的噪声
ar = [1, 1/2]; % AR模型
v_ar = filter(1, ar, v);
% 计算噪声的自相关函数
[acf, lags] = xcorr(v_ar, order_ma, 'unbiased'); % 计算自相关函数

% 只取正向自相关
acf = acf(order_ma+1:end); 

% 设计一个矩阵，使用最小二乘法估算MA系数
R = toeplitz(acf(1:end-1));  % 自相关矩阵（长度为阶数-1）
b = acf(2:end);               % 自相关向量（去除自相关函数的第一个元素）

% 使用最小二乘法求解MA系数
ma_coeff = R\b;

% 打印结果
disp('估算的MA系数：');
disp(ma_coeff);

% 使用估算的MA系数生成去噪后的信号
v_ma_filtered = filter([1; -ma_coeff], 1, v_ar); 

% 绘制结果
figure;
subplot(3,1,1);
plot(v); 
title('原始噪声信号');
xlabel('时间');
ylabel('幅度');
grid on;

subplot(3,1,2);
plot(v_ma_filtered); 
title('去噪后的信号（MA滤波）');
xlabel('时间');
ylabel('幅度');
grid on;

subplot(3,1,3);
stem(ma_coeff); 
title('估算的MA系数');
xlabel('系数编号');
ylabel('系数值');
grid on;
