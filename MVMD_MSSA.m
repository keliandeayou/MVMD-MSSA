clc
clear
close all
% 读取数据
load('SJ.mat');
fs=6000;
% 读取前10000长度的信号
len=6000;
%DATAt=SJ031(1:6000,1);
s22=SJ031(1:6000,2);
s33=SJ031(1:6000,3);
s2=SJ031(1:6000,2:3);
% 采样时间
t = (0:len-1)/fs;
alpha=1107;
K=8;
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度
%MVMD分解
[imf2, ~, ~] = MVMD(s2, alpha, tau, K, DC, init, tol);    
% imf2 =imf2';
% u1 = imf2(:, 1); % 分解后的第一个分量
% u2 = imf2(:, 2); % 分解后的第二个分量
% u3 = imf2(:, 3); % 分解后的第三个分量
% u4 = imf2(:, 4); % 分解后的第四个分量
% u5 = imf2(:, 5); % 分解后的第五个分量
% u6 = imf2(:, 6); % 分解后的第六个分量
% u7 = imf2(:, 7); % 分解后的第七个分量
% u8 = imf2(:, 8); % 分解后的第七个分量
imf2 = permute(imf2, [2, 1, 3]); 
column2_imf = imf2(:, :, 1);
column3_imf = imf2(:, :, 2);
%向量2.3.5相关系数较大
%%计算合成数据的傅里叶谱
% 假设 u1, u2, ..., u6 分别代表分解后的6个IMF分量
u1 = column2_imf(:, 1); % 分解后的第一个分量
u2 = column2_imf(:, 2); % 分解后的第二个分量
u3 = column2_imf(:, 3); % 分解后的第三个分量
u4 = column2_imf(:, 4); % 分解后的第四个分量
u5 = column2_imf(:, 5); % 分解后的第五个分量
u6 = column2_imf(:, 6); % 分解后的第六个分量
u7 = column2_imf(:, 7); % 分解后的第七个分量
u8 = column2_imf(:, 8); % 分解后的第七个分量

u11 = column3_imf(:, 1); % 分解后的第一个分量
u22 = column3_imf(:, 2); % 分解后的第二个分量
u33 = column3_imf(:, 3); % 分解后的第三个分量
u44 = column3_imf(:, 4); % 分解后的第四个分量
u55 = column3_imf(:, 5); % 分解后的第五个分量
u66 = column3_imf(:, 6); % 分解后的第六个分量
u77 = column3_imf(:, 7); % 分解后的第七个分量
u88 = column3_imf(:, 8); % 分解后的第七个分量
imf5=[u5,u55];
imf6=[u6,u66];

% 设置采样频率和采样点数
fs = 6000; % 采样率
N = size(u1, 1); % 采样点数
% 计算频率轴
f = (-N/2:N/2-1) * fs / N; % 频率轴


%对IMF1进行MSSA分解
DATA = u1; % 输入数据
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析Hankel矩阵的秩数
flow =lower_freq1;         %分量最低频率
fhigh=upper_freq1;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f1,~] = mssa_2d(DATA,dt,P,flow,fhigh,meth);

%%画出滤波后的频谱图

% % 计算频谱
% Fs = 6000;  % 采样频率
% L = size(DATA_f1, 1);  % 数据长度
% f = Fs*(0:(L/2))/L;  % 频率范围
% 
% % 对经过MSSA滤波后的数据应用频谱分析
% DATA_f_fft = fft(DATA_f1);
% P2 = abs(DATA_f_fft/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% % 绘制频谱图
% figure();
% plot(f, P1);
% title('经过MSSA滤波后的频谱图');
% xlabel('频率 (Hz)');
% ylabel('能量');

%对IMF2进行MSSA分解
DATA = u2; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析的总时间窗长度
flow = lower_freq2;         %分量最低频率
fhigh=upper_freq2;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f2,~] = mssa_2d(DATA,dt,P,flow,fhigh,meth);


%对IMF3进行MSSA分解
DATA = u3; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P =1;           % MSSA分析的总时间窗长度
flow = lower_freq3;         %分量最低频率
fhigh=upper_freq3;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f3,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);



%对IMF4进行MSSA分解
DATA = u4; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P =1;           % MSSA分析的总时间窗长度
flow = lower_freq4;         %分量最低频率
fhigh=upper_freq4;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f4,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);



%对IMF5进行MSSA分解
DATA5 = imf5;        
dt = 0.17 / 1000;   % 数据采样时间间隔
P = 1;              % hankel矩阵的秩
flow=lower_freq5;         %分量最低频率
fhigh=upper_freq5;            % 分量最高频率
meth=1;             % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f5,sing] = mssa_2d(DATA5,dt,P,flow,fhigh,meth);
DATA_f15 = DATA_f5(:, 1); % 分解后的第一个分量
DATA_f25 = DATA_f5(:, 2); % 分解后的第二个分量
%%画出滤波后的频谱图

%对IMF6进行MSSA分解
DATA6 = imf6; % 噪声容忍度
dt = 0.17 / 1000;    % 数据采样时间间隔
P = 1;               % MSSA分析的总时间窗长度
flow=lower_freq5;         %分量最低频率
fhigh=upper_freq5;            % 分量最高频率
meth=1;              % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f6,sing] = mssa_2d(DATA6,dt,P,flow,fhigh,meth);
DATA_f16 = DATA_f6(:, 1); % 分解后的第一个分量
DATA_f26 = DATA_f6(:, 2); % 分解后的第二个分量


%对IMF7进行MSSA分解
DATA = u7; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析的总时间窗长度
flow=lower_freq7;         %分量最低频率
fhigh=upper_freq7;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f7,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);

%对IMF8进行MSSA分解
DATA = u8; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析的总时间窗长度
flow=lower_freq8;         %分量最低频率
fhigh=upper_freq8;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f8,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);



%%信号重构
reconstructed_signal2=u1+u22+u3+u4+u5+u6+u7+u8;
reconstructed_signal22=u5+u6;
reconstructed_signal222=DATA_f15+DATA_f16;
%reconstructed_signal222=DATA_f15+DATA_f16;

reconstructed_signal3=u11+u22+u33+u44+u55+u66+u77+u88;
reconstructed_signal33=u55+u66;
reconstructed_signal333=DATA_f25+DATA_f26;
%reconstructed_signal333=DATA_f25+DATA_f26;

%%互相关时延
x=s2;
y=s3;
x=reconstructed_signal22;
y=reconstructed_signal33;
Fs=6000;
corr = xcorr(x, y); 
%corr = xcorr(x, y, 'coeff'); 
M = length(corr);
dt = (-floor(M/2):floor(M/2))/Fs;
[max_corr, max_pos] = max(abs(corr));
delay_time = dt(max_pos);
fprintf('时间延迟为 %.2f 秒（或者 %d 个样本）\n', delay_time, round(delay_time*Fs));
dt =dt';
% 绘制互相关图
figure;
plot(dt, corr);
xlabel('延迟（秒）');
ylabel('互相关函数');
title('两个重构信号的互相关图');
c = 1478; % m/s
% 根据时间延迟计算泄漏点距离传感器的位置
d=6.6; % m
d2 = (d+c * delay_time) / 2;
