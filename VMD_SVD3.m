clc
clear
close all
% 读取数据
load('SJ.mat');

fs=6000;
% 读取前10000长度的信号
len=6000;
s=SJ0652(1:6000,3); 
% s3=SJ014(1:6000,3);
% 采样时间
t = (0:len-1)/fs;
alpha=2000;
K=8;
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度
% data2=s2;             %data为待分解数据，可设置为带分解数据集中的任意一条，根据需求可以直接编写一个循环实现所有数据的分解
% data3=s3;             %data为待分解数据，可设置为带分解数据集中的任意一条，根据需求可以直接编写一个循环实现所有数据的分解


[u, ~, ~] = VMD(s, alpha, tau, K, DC, init, tol);%采用最优参数进行分解
% [u, ~, ~] = MVMD(s, alpha, tau, K, DC, init, tol);

v2 = u';             % 转置得到列向量
% 分别用u1、u2、u3......u7来表示这7个分量
u1 = v2(:, 1);
u2 = v2(:, 2);
u3 = v2(:, 3);
u4 = v2(:, 4);
u5 = v2(:, 5);
u6 = v2(:, 6);
u7 = v2(:, 7);
u8 = v2(:, 8);
reconstructed_signal3=u1+u2+u3+u4+u5+u6+u7+u8;
reconstructed_signal33=u5+u6;
% 假设已将数据分解为7个IMF分量u1、u2、u3......u7
% 计算每个分量的峭度值
kurt1 = kurtosis(u1);
kurt2 = kurtosis(u2);
kurt3 = kurtosis(u3);
kurt4 = kurtosis(u4);
kurt5 = kurtosis(u5);
kurt6 = kurtosis(u6);
kurt7 = kurtosis(u7);

% 找出峭度值较大的前2个分量
% [~, idx] = sort([kurt1, kurt2, kurt3, kurt4, kurt5, kurt6], 'descend');
[~, idx] = sort([kurt1, kurt2, kurt3, kurt4, kurt5, kurt6, kurt7], 'descend');
imf_select = idx(1:2);

% 对选定的IMF分量进行重构
v2_select = v2(:, imf_select);
s_recon = sum(v2_select, 2);

% 绘制重构后的信号
figure;
plot(t, s_recon);
xlabel('Time');
ylabel('Amplitude');
title('Reconstructed signal with selected IMF components');


% 计算重构信号的Hankel矩阵
L = length(s_recon);
K1 = 2510; % 选取的Hankel矩阵的列数
H = hankel(s_recon(1:K1),s_recon(K1:L));

% 对Hankel矩阵进行奇异值分解
[U, S, V] = svd(H);

% 计算奇异值的差分谱
diffS = abs(diff(diag(S)));

% 单边极大值原则选择差分谱峰值，确定有效秩阶次
n = find(diffS == max(diffS), 1, 'first');
r = n + 1;

% 根据选择的秩进行信号重构
U_r = U(:, 1:r);
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
s_recon_filter = U_r * S_r * V_r';

% 计算阈值
threshold = mean(H(:)); % 或者使用中值，如：median(H(:));

% 进行信号滤波
s_filter = medfilt1(s_recon_filter, 3); % 或者使用其他滤波器，如：sgolayfilt(s_recon_filter, 3, 5);
s_filter(s_filter<threshold) = 0;

s_filter_transpose = s_filter';
s_filter_col = s_filter_transpose(:,1);

s_recon_col = s_recon(:,1);
s_recon_col = s_recon_col(1:3491);


% 计算信噪比
noise = s_recon_col- s_filter_col;
SNR = 10*log10(norm(s_recon)^2/norm(noise)^2);

% 重构信号的时域特性图
figure;
plot(s_recon_col);
xlabel('Sample');
ylabel('Amplitude');
title('Time Domain Characteristics of Reconstructed Signal');

% 重构信号的频域特性图
figure;
Fs = 6000; % 采样率, 假设为6000Hz
L = length(s_recon_col); % 信号长度
Y = fft(s_recon_col); % 进行FFT变换
P2 = abs(Y/L); % 计算频率幅值（除以信号长度）
P1 = P2(1:L/2+1); % 取正半轴
P1(2:end-1) = 2*P1(2:end-1); % 对于奇数长度的信号，需要对中间值进行倍增
f = Fs*(0:(L/2))/L; % 计算频率轴
plot(f,P1);
title('Single-Sided Amplitude Spectrum of Reconstructed Signal');
xlabel('f (Hz)');
ylabel('|P1(f)|');