clc
clear
close all
% 读取数据
load('SJ.mat');

fs=6000;
% 读取前10000长度的信号
len=6000;
s2=SJ072(1:6000,3);
% s3=SJ072(1:6000,3);
% 采样时间
t = (0:len-1)/fs;
alpha=2000;
K=6;
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度
 
[u2, ~, ~] = MVMD(s2, alpha, tau, K, DC, init, tol);
u2 = u2'; 
% IMF分量相加重构信号
reconstructed_signal = sum(u2,2);


% 计算原始信号和重构信号的平方和
signal_power = sum(s2.^2);
error_power = sum((s2 - reconstructed_signal).^2);

% 计算SNR0
SNR0 = 10 * log10(signal_power / error_power);