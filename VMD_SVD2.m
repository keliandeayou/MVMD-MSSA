%% 改进的基于麻雀搜索算法优化变分模态分解参数
clc
clear
close all
% 读取数据
load('SJ.mat');

fs=6000;
% 读取前10000长度的信号
len=6000;
s=SJ0551(1:6000,2);
% s3=SJ031(1:6000,3);

% 采样时间
t = (0:len-1)/fs;0;


%% 设定改进的麻雀搜索算法参数
popsize =20;   % 种群大小,可更改
iter = 100;   % 最大迭代次数,可更改
dim = 2;   % VMD变量个数
lb = [800 3];   % alpha范围 K范围   下限
ub = [2000 8];  % 上限
%% 改进的麻雀搜索算法SSA优化VMD参数
tic ,  % 开始计时
trainlabel=1;%是否重新训练。1就是重新训练，0就是直接调取上次寻优得到的结果
if trainlabel==1
    [~,alpha,K]=IAMSSA_VMD(popsize,iter,dim,lb,ub,1);  % 0表示不保存IMF，1表示导出IMF并保存
    save result alpha K 
else
    load result alpha K 
end
toc,  % 结束计时



clc
clear
close all
% 读取数据
load('SJ.mat');

fs=6000;
% 读取前10000长度的信号
len=6000;
s=SJ0652(1:6000,2);
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

%[u, ~, ~] = MVMD(s, alpha, tau, K, DC, init, tol);

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
reconstructed_signal2=u1+u2+u3+u4+u5+u6+u7+u8;
reconstructed_signal22=u5+u6;

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
% [~, idx] = sort([kurt1, kurt2, kurt3, kurt4, kurt5, kurt6 ], 'descend');
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
K1 = 2510;            % 选取的Hankel矩阵的列数
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

%对信号进行规范化调整
s_filter_transpose = s_filter';
s_filter_col = s_filter_transpose(:,1);
% s_filter_col = s_filter_col(1:3001);

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





%计算各分量与原信号相关系数
v2 = u2';             % 转置得到列向量
v3 = u3';             % 转置得到列向量
kr2 = zeros(1, 7);   % 用于存储相关系数的向量
for i = 1:7
    [Ra2, ~] = xcorr(v2(:, i), s2);
    [Rb2, lag2] = xcorr(s2);% 计算自相关函数Rb
    C2 = corrcoef(Ra2, Rb2);
    kr2(i) = C2(1, 2);
end
% 画出柱状图
bar(kr2);
title('各分量与原信号相关系数');
xlabel('分量');
ylabel('相关系数');

kr3 = zeros(1, 7);   % 用于存储相关系数的向量
for i = 1:7
    [Ra3, ~] = xcorr(v3(:, i), s3);
    [Rb3, lag3] = xcorr(s3);% 计算自相关函数Rb
    C3 = corrcoef(Ra3, Rb3);
    kr3(i) = C3(1, 2);
end
% 画出柱状图
bar(kr3);
title('各分量与原信号相关系数');
xlabel('分量');
ylabel('相关系数');




x=imf_list_col;
y=imf_list_col3;
Fs=6000;
% 假设两个信号数据分别为 s1 和 s2，采样频率为 Fs，数据长度分别为 N1 和 N2
% 需要提前确保 N1 和 N2 大小一致
% 计算信号 s1 和 s2 的互相关函数
corr = xcorr(x, y);
% 获取互相关函数的长度
M = length(corr);
% 计算时间延迟（单位：秒）
% 对于 Fs = 1，时间延迟的单位是秒；对于 Fs 不等于 1，时间延迟的单位是样本数
dt = (-floor(M/2):floor(M/2))/Fs;
% 找出互相关函数的最大值及其位置，即为两个信号之间的最大相关值和时间延迟
[max_corr, max_pos] = max(abs(corr));
delay_time = dt(max_pos);
% 输出结果
fprintf('时间延迟为 %.2f 秒（或者 %d 个样本）\n', delay_time, round(delay_time*Fs));

% 绘制互相关图
figure;
plot(lag, coeff);
title('X和Y的互相关');
xlabel('时间/s'); 
ylabel('相关系数');

% 管道泄漏的声信号沿管道的传播速度
c = 1478; % m/s
% 根据时间延迟计算泄漏点距离传感器的位置
d=6.6; % m
d2 = (d+c * delay_time) / 2;

