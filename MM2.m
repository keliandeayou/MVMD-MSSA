clc
clear
close all
% 读取数据
load('SJ.mat');
fs=6000;
% 读取前10000长度的信号
len=6000;
% s3=SJ072(1:6000,2);
s3=SJ041(1:6000,3);
% 采样时间
t = (0:len-1)/fs;
alpha=1307;
K=8;
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度
% [imf2, ~, ~] = VMD(s3, alpha, tau, K, DC, init, tol);
%MVMD分解
[imf2, ~, ~] = MVMD(s3, alpha, tau, K, DC, init, tol);
imf2 =imf2'; 
%%计算合成数据的傅里叶谱
% 假设 u1, u2, ..., u6 分别代表分解后的6个IMF分量
u1 = imf2(:, 1); % 分解后的第一个分量
u2 = imf2(:, 2); % 分解后的第二个分量
u3 = imf2(:, 3); % 分解后的第三个分量
u4 = imf2(:, 4); % 分解后的第四个分量
u5 = imf2(:, 5); % 分解后的第五个分量
u6 = imf2(:, 6); % 分解后的第六个分量
u7 = imf2(:, 7); % 分解后的第七个分量
u8 = imf2(:, 8); % 分解后的第七个分量
% u1 = u11; % 分解后的第一个分量
% u2 = u22; % 分解后的第二个分量
% u3 = u33; % 分解后的第三个分量
% u4 = u44; % 分解后的第四个分量
% u5 = u55; % 分解后的第五个分量
% u6 = u66; % 分解后的第六个分量
% u7 = u77; % 分解后的第七个分量
% u8 = u88; % 分解后的第七个分量
% 设置采样频率和采样点数
fs = 6000; % 采样率
N = size(u1, 1); % 采样点数

% 计算频率轴
f = (-N/2:N/2-1) * fs / N; % 频率轴
chonggou=u5+u6;
% 计算峭度值
kurtosis_u1 = kurtosis(u1);
kurtosis_u2 = kurtosis(u2);
kurtosis_u3 = kurtosis(u3);
kurtosis_u4 = kurtosis(u4);
kurtosis_u5 = kurtosis(u5);
kurtosis_u6 = kurtosis(u6);
kurtosis_u7 = kurtosis(u7);
kurtosis_u8 = kurtosis(u8);




%对IMF1进行MSSA分解
DATA = u1; % 输入数据
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析的总时间窗长度（Hankel矩阵的秩数）
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
P = 1;           % MSSA分析的总时间窗长度
flow = lower_freq3;         %分量最低频率
fhigh=upper_freq3;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f3,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);



%对IMF4进行MSSA分解
DATA = u4; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析的总时间窗长度
flow = lower_freq4;         %分量最低频率
fhigh=upper_freq4;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f4,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);



%对IMF5进行MSSA分解
DATA = u5; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析的总时间窗长度
flow=lower_freq5;         %分量最低频率
fhigh=upper_freq5;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f5,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);

%%画出滤波后的频谱图




%对IMF6进行MSSA分解
DATA = u6; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析的总时间窗长度
flow=lower_freq6;         %分量最低频率
fhigh=upper_freq6;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f6,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);

%%画出滤波后的频谱图



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

DATA = u8; % 噪声容忍度
dt = 0.17 / 1000;     % 数据采样时间间隔
P = 1;           % MSSA分析的总时间窗长度
flow=lower_freq7;         %分量最低频率
fhigh=upper_freq7;            % 分量最高频率
meth=1;            % 收敛准则容忍度
%DATA_f: 执行MSSA分析后的各个频率分量（特征模态或轨道）组成的二维矩阵。每一行对应一个频率分量，每一列对应一个时间步骤。
%sing: MSSA的奇异值谱图，包含信号的频率和能量信息。
[DATA_f8,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);


%%信号重构
reconstructed_signal3=u1+u2+u3+u4+u5+u6+u7+u8;
reconstructed_signal33=u5+u6;
reconstructed_signal333=DATA_f5+DATA_f6;

%将滤波后的 IMF 分量相加，并进行逆变换，你可以得到恢复后的有效泄露数据。
recovered_data2 = ifft(reconstructed_signal2);
recovered_data = ifft(reconstructed_signal);

