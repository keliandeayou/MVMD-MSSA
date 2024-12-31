%% 改进的基于麻雀搜索算法优化变分模态分解参数
    clc
    clear
    close all
    % 读取数据
    load('SJ.mat');
    
    fs=6000;
    % 读取前10000长度的信号
    len=6000;
    s=SJ041(1:6000,2);
    % s3=SJ031(1:6000,3);
    
    % 采样时间
    t = (0:len-1)/fs;0;
    
    
    %% 设定改进的麻雀搜索算法参数
    popsize =25;   % 种群大小,可更改
    iter = 25;   % 最大迭代次数,可更改
    dim = 2;   % VMD变量个数
    lb = [1000 6];   % alpha范围 K范围下限
    ub = [3000 10];  % 上限
    % %麻雀优化算法
    % [fMin , bestX,Convergence_curve1 ] = SSA(popsize, iter,lb,ub,dim,@(x)fun(x,s));
    % %基础粒子群优化算法
    % [gBestScore,gBest,cg_curve]=PSO(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % %狼
    % [Leader_score,Leader_pos,Convergence_curve]=WOA(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % %WSO
    % [fmin0,gbest,ccurve]=WSO(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % %MFO
    % [Best_flame_score,Best_flame_pos,Convergence_curve]=MFO(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % %
    % [Alpha_score,Alpha_pos,Convergence_curve]=GWO(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % flag=0;   % flag为1，则添加文本水印，为0，不添加水印
    % fig_show(1)
    % plot(Convergence_curve1,'k-.o','linewidth',1)
    % legend('最佳适应度')
    % grid on
    % xlabel('进化代数')
    % ylabel('包络熵')
    % title('IAMSSA进化曲线')
    % fig_text(flag)
    % 
    % %导入EXCEL
    % writetable(Convergence_curve, "1.csv")
    %% 改进的麻雀搜索算法SSA优化VMD参数
    tic ,  % 开始计时
    trainlabel=1;%是否重新训练。1就是重新训练，0就是直接调取上次寻优得到的结果
    if trainlabel==1
        [~,alpha,K]=IAMSSA_MVMD(popsize,iter,dim,lb,ub,1);  % 0表示不保存IMF，1表示导出IMF并保存
        save result alpha K 
    else
        load result alpha K 
    end
toc,  % 结束计时
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度
[u, ~, ~] = VMD(s, 1109.1, tau, 8, DC, init, tol);

fs=6000;
% 读取前10000长度的信号
len=6000;
s2=SJ072(1:6000,2);
s3=SJ072(1:6000,3); 
% s3=SJ014(1:6000,3);
% 采样时间
t = (0:len-1)/fs;

%互谱分析
s2_pre = wiener2(s2); % 对信号s2进行去噪处理
s3_pre = wiener2(s3); % 对信号s3进行去噪处理

fs = 6000; % 采样频率，单位Hz
N = length(s2_pre);
S2_fft = fft(s2_pre, N)/N; % 使用MATLAB中的fft函数进行傅里叶变换
S3_fft = fft(s3_pre, N)/N; % 使用MATLAB中的fft函数进行傅里叶变换

f = fs*(0:(N/2))/N; % 计算频率向量，单位Hz
S2_psd = 2*abs(S2_fft(1:N/2+1)).^2; % 计算信号的功率谱密度
S3_psd = 2*abs(S3_fft(1:N/2+1)).^2; % 计算信号的功率谱密度

[P_S2S3, f] = cpsd(s2_pre, s3_pre, [], [], N, fs); % 计算信号的互谱密度
P_S2S3_norm = abs(P_S2S3)./sqrt(S2_psd.*S3_psd); % 归一化处理，得到互谱的相关程度

omega_s1 = 1000;  % rad/s
omega_s2 = 1000;  % rad/s
omega = (2*pi*f)';
index = find((omega>=omega_s1)&(omega<=omega_s2)); % 查找估计特征频带 ωs1～ωs2 的索引
P_S2S3_norm_filt = P_S2S3_norm(index); % 从归一化后的互谱中提取估计特征频带 ωs1～ωs2 的信息



%% 采用优化所得参数进行VMD分解
% load result alpha K %加载最优参数
% alpha=2000;
% K=7;
clc
clear
close all
% 读取数据
load('SJ.mat');

fs=6000;
% 读取前10000长度的信号 
len=6000;
s2=SJ072(1:6000,2);
s3=SJ072(1:6000,3); 

% 采样时间
t = (0:len-1)/fs;
% alpha=2000;
% K=7;
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度
data2=s2;             %data为待分解数据，可设置为带分解数据集中的任意一条，根据需求可以直接编写一个循环实现所有数据的分解
data3=s3;             %data为待分解数据，可设置为带分解数据集中的任意一条，根据需求可以直接编写一个循环实现所有数据的分解

% [u2, ~, ~] = MVMD(s2, alpha, tau, K, DC, init, tol);
% [u3, ~, ~] = MVMD(s3, alpha, tau, K, DC, init, tol);
[u2, ~, ~] = VMD(s2, 605, tau, 3, DC, init, tol);%采用最优参数进行分解
[u3, ~, ~] = VMD(s3, 500, tau, 3, DC, init, tol);%采用最优参数进行分解
%输出变量 u是一个矩阵，其中每一行代表一个IMF分量，每一列代表一个时间点。
%K1+1是IMF的分解个数，实际选取分解的层数为K1，也就是VMD分解得到的IMF分量中有K1个分量需要保留，其余的分量舍去，即选取前K1个IMF分量进行重构，s2_vmd = sum(IMF(1:K1, :))。


%计算各分量与原信号相关系数
v2 = u2';             % 转置得到列向量
v3 = u3';             % 转置得到列向量
kr2 = zeros(1, 3);   % 用于存储相关系数的向量
for i = 1:3
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

kr3 = zeros(1, 3);   % 用于存储相关系数的向量
for i = 1:3
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


% 选取相关系数大于最大相关系数一半的分量进行重构
max_corr2 = max(kr2);
threshold2 = max_corr2/3;
xr2 = zeros(size(v2,1),1);

for i = 1:size(v2,2)
    if kr2(i) > threshold2
        xr2 = xr2 + v2(:,i);
    end
end

% 选取相关系数大于最大相关系数一半的分量进行重构
max_corr3 = max(kr3);
threshold3 = max_corr3/3;
xr3 = zeros(size(v3,1),1);

for i = 1:size(v3,2)
    if kr3(i) > threshold3
        xr3 = xr3 + v3(:,i);
    end
end


% 输出重构后的向量矩阵
disp('重构后的向量矩阵：');
disp(xr2');

% 输出重构后的向量矩阵
disp('重构后的向量矩阵：');
disp(xr3');


%重构信号与原始信号对比图
t2 = 1:length(s2); %生成时间轴
figure;
plot(t2,s2,'b');
hold on;
plot(t2,xr2,'r');
title('原始信号和重构信号对比');
xlabel('时间');
ylabel('幅值');
legend('原始信号','重构信号');
%两个对比图分开的
t2 = 1:length(s2); %生成时间轴
figure;
subplot(2,1,1);
plot(t2,s2);
title('原始信号');
xlabel('时间');
ylabel('幅值');
subplot(2,1,2);
plot(t2,xr2);
title('重构信号');
xlabel('时间');
ylabel('幅值');



%重构信号与原始信号对比图
t3 = 1:length(s3); %生成时间轴
figure;
plot(t3,s3,'b');
hold on;
plot(t3,xr3,'r');
title('原始信号和重构信号对比');
xlabel('时间');
ylabel('幅值');
legend('原始信号','重构信号');
%两个对比图分开的
t3 = 1:length(s3); %生成时间轴
figure;
subplot(2,1,1);
plot(t3,s3);
title('原始信号');
xlabel('时间');
ylabel('幅值');
subplot(2,1,2);
plot(t3,xr3);
title('重构信号');
xlabel('时间');
ylabel('幅值');


x=s_filter_col;
y=s_filter_col3;
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


% %互相关分析
% x=xr2;
% y=xr3;
% [coeff, lag] = xcorr(x, y);
% [~, max_index] = max(coeff);
% delay = lag(max_index);
% sampling_rate = 6000; % 假设采样率为 6000
% delay_time = delay / sampling_rate;


% % 绘制互相关图
% figure;
% plot(lag, coeff);
% title('X和Y的互相关');
% xlabel('时间/s');
% ylabel('相关系数');


% 管道泄漏的声信号沿管道的传播速度
c = 1478; % m/s
% 根据时间延迟计算泄漏点距离传感器的位置
d=6.6; % m
d2 = (d+c * delay_time) / 2;







 

 
 
 %% 假设u为分量矩阵，其中N为分量的数目，而M为每个分量的长度
%% 初始化相关系数矩阵，每个分量之间的相关系数按照矩阵对称性对称取值
corr_mat = corrcoef(u);
% 输出相关系数矩阵
disp('分量之间的相关系数矩阵：');
disp(corr_mat);


% 绘制相关性矩阵
figure;
imagesc(abs(corr_mat));
colorbar;
colormap(jet);
title('分量之间的相关系数');
xlabel('分量编号');
ylabel('分量编号');





kurt_vals = kurtosis(u');

kurt_thresh = 2; % 假设选择峭度值大于2的模态分量进行重构
important_modes = find(kurt_vals > kurt_thresh);
reconstructed_signal = sum(u(:, important_modes), 2);

plot(x);
hold on;
plot(reconstructed_signal);
xlabel('Time');
ylabel('Amplitude');
legend('Original Signal', 'Reconstructed Signal');



% average=mean(omega);
% plot(omega);
% 
% 
% figure('Name','频谱图','Color','white');
% for i = 1:K
% p=abs(fft(u(:,i))); %并fft，得到p，就是包络线的fft---包络谱 
% subplot(K,1,i);
% plot((0:len-1)*fs/len,p) %绘制包络谱
% xlim([0 fs/2]) %展示包络谱低频段，这句代码可以自己根据情况选择是否注释
% if i ==1
% title('频谱图'); xlabel('频率'); ylabel(['IMF' int2str(i)]);%int2str(i)是将数值i四舍五入后转变成字符，y轴命名
% else
% xlabel('频率'); ylabel(['IMF' int2str(i)]);%int2str(i)是将数值i四舍五入后转变成字符，y轴命名
% end
% end
% set(gcf,'color','w');
%% 作图
% 画imf分量的时域图
[m,~]=size(u);
m=m+1;
figure
for i=1:m-1
    subplot(m,1,i)
    plot(t,u(i,:),'b-','linewidth',1);hold on
    ylabel(['IMF',num2str(i)]);
end
res = s'-sum(u,1);hold on
subplot(m,1,m)
plot(t,res,'b-','linewidth',1)
ylabel('Res');
xlabel('t/s')
hold off
title('VMD分解后的时域图')
% 画imf分量的频域图
figure
for i=1:m-1
    subplot(m,1,i)
    %% FFT 变换
    [cc,y_f]=plot_fft(u(i,:),fs,1);%对输入的时域信号u(i,:)进行FFT变换
    plot(y_f,cc,'b','LineWIdth',1);hold on
    ylabel(['FFT of IMF',num2str(i)]);
end
hold on
subplot(m,1,m)
[cc,y_f]=plot_fft(res,fs,1);
plot(y_f,cc,'b','LineWIdth',1);hold off
ylabel('FFT of Res');
xlabel('f/Hz')
fig_text(flag)
title('VMD分解后的频域图')
%% 未优化VMD分解结果
alpha=1000;  % 惩罚因子，也称平衡参数
K=5;  % 分解的模态数
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度
%% 未优化VMD分解
% u：分解模式的集合
% u_hat：模式的频谱
% omega：估计模式中心频率
[u, u_hat, omega] = VMD(s, alpha, tau, K, DC, init, tol);
%% 作图
% 画imf分量的时域图
[m,~]=size(u);
m=m+1;
figure
for i=1:m-1
    subplot(m,1,i)
    plot(t,u(i,:),'b-','linewidth',1);hold on
    ylabel(['IMF',num2str(i)]);
end
res = s'-sum(u,1);hold on
subplot(m,1,m)
plot(t,res,'b-','linewidth',1)
ylabel('Res');
xlabel('t/s')
hold off
title('VMD分解后的时域图')
% 画imf分量的频域图
figure
for i=1:m-1
    subplot(m,1,i)
    %% FFT 变换
    [cc,y_f]=plot_fft(u_hat(i,:),fs,1);
    plot(y_f,cc,'b','LineWIdth',1);hold on
    ylabel(['FFT of IMF',num2str(i)]);
end
hold on
subplot(m,1,m)
[cc,y_f]=plot_fft(res,fs,1);
plot(y_f,cc,'b','LineWIdth',1);hold off
ylabel('FFT of Res');
xlabel('f/Hz')
fig_text(flag)
title('VMD分解后的频域图')
tic
clc
clear all




clc
clear
close all
% 读取数据
load('SJ.mat');

fs=6000;
% 读取前6000长度的信号
len=6000;
s=SJ110(1:6000,2);
% 采样时间
t = (0:len-1)/fs;
% 调用MVMD函数进行分解
[numIMFs,IMF] = MVMD(s);
%2.构建MVMD分解器
vmd = matvmd(s);

%3.进行MVMD分解
[Coef,v] = vmd.vmd();

%4.选择合适的分量
energy = sum(abs(Coef).^2,2);
energy_ratio = energy./sum(energy);
components = v(:,energy_ratio>0.1);

%绘制分解结果
figure;
subplot(2,1,1);
plot(signal);
title('Original Signal');
subplot(2,1,2);
plot(sum(Coef,1)');
title('Sum of Coefficients');

