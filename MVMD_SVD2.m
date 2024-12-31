clc
clear
close all
% 读取数据
load('SJ.mat');

fs=6000;
% 读取前10000长度的信号
len=6000;
s2=SJ072(1:6000,2);
% s3=SJ072(1:6000,3);
% 采样时间
t = (0:len-1)/fs;
alpha=1597;
K=7;
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度
% data2=s2;             %data为待分解数据，可设置为带分解数据集中的任意一条，根据需求可以直接编写一个循环实现所有数据的分解
% data3=s3;             %data为待分解数据，可设置为带分解数据集中的任意一条，根据需求可以直接编写一个循环实现所有数据的分解
% 

% 
%  [u2, ~, ~] = MVMD(s2, alpha, tau, K, DC, init, tol);

% 时域分析
% 计算信号的均值和标准差，确定泄漏信号的时域特征
s2_mean = mean(s2);%平均值
s2_std = std(s2);%标准差
% s3_mean = mean(s3);
% s3_std = std(s3);
fprintf('s2(t)的平均值: %.2f, 标准差: %.2f\n', s2_mean, s2_std);
% fprintf('s3(t)的平均值: %.2f, 标准差: %.2f\n', s3_mean, s3_std);

% 频域分析
% 对信号进行FFT变换
L = length(s2);
f = fs*(0:(L/2))/L;

%对s2频域分析
Y2 = fft(s2);%计算信号s2在频域的傅里叶变换，得到频域复信号Y2。
P2_2 = abs(Y2/L);%计算傅里叶变换后频率响应的振幅。
P1_2 = P2_2(1:L/2+1);%表示取出 P2_2 的前一半，并将其存储在P1_2 中。
P1_2(2:end-1) = 2*P1_2(2:end-1); %通过将频率响应乘以 2，将复信号的幅值翻倍。
% % 画出波形谱图
% figure();
% plot(f,P1_2) 
% title('Single-Sided Amplitude Spectrum of s2(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')


%对s3频域分析
% Y3 = fft(s3);%计算信号s3在频域的傅里叶变换，得到频域复信号Y2。
% P2_3 = abs(Y3/L);%计算傅里叶变换后频率响应的振幅。
% P1_3 = P2_3(1:L/2+1);%表示取出 P2_3 的前一半，并将其存储在P1_3 中。
% P1_3(2:end-1) = 2*P1_3(2:end-1); %通过将频率响应乘以 2，将复信号的幅值翻倍。
% % 画出波形谱图
% figure();
% plot(f,P1_3) 
% title('Single-Sided Amplitude Spectrum of s3(t)')
% xlabel('f (Hz)')
% ylabel('|P2(f)|')



% %对无泄漏信号进行分析
% x1=bulou1(1:6000,2);
% x2=bulou1(1:6000,3);
% % 时域分析
% % 比较信号的均值和标准差等，确定泄漏信号的时域特征
% mean_x1 = mean(x1);
% std_x1 = std(x1);
% mean_x2 = mean(x2);
% std_x2 = std(x2);
% fprintf('x1(t)的平均值: %.2f, 标准差: %.2f\n', mean_x1, std_x1);
% fprintf('x2(t)的平均值: %.2f, 标准差: %.2f\n', mean_x2, std_x2);
% 
% % 频域分析
% % 对信号进行FFT变换
% L = length(x1);
% Y1 = fft(x1);
% Y2 = fft(x2);
% P2_1 = abs(Y1/L);
% P1_1 = P2_1(1:L/2+1);
% P1_1(2:end-1) = 2*P1_1(2:end-1); 
% P2_2 = abs(Y2/L);
% P1_2 = P2_2(1:L/2+1);
% P1_2(2:end-1) = 2*P1_2(2:end-1); 
% f = fs*(0:(L/2))/L;
% 
% % 画出波形谱图
% figure();
% plot(f,P1_1) 
% title('Single-Sided Amplitude Spectrum of X1(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% figure();
% plot(f,P1_2) 
% title('Single-Sided Amplitude Spectrum of X2(t)')
% xlabel('f (Hz)')
% ylabel('|P2(f)|')


[u2, ~, ~] = MVMD(s2, alpha, tau, K, DC, init, tol);

     
% u2 = double(u2);

% 定义参数
lower = 1600;
upper = 2000;

L2=6000;
% L3=6000;


% 计算特征频带内的泄漏信号（s2）
P1_2 = abs(Y2/L2).^2;
freq = (0:length(P1_2)-1)*fs/L2;
chosen_freq = (freq >= lower) & (freq <= upper);
P1_2_s1 = P1_2(chosen_freq);
% % 计算特征频带内的泄漏信号（s3）
% P1_3 = abs(Y3/L3).^2;
% freq = (0:length(P1_3)-1)*fs/L3;
% 
% chosen_freq = (freq >= lower) & (freq <= upper);
% P1_3_s1 = P1_3(chosen_freq);


%计算各分量与原信号相关系数
u2 = u2'; 
rho_list = zeros(1, 7); % 创建一个长度为7的零向量
for i = 1:7 % 假设分量矩阵u2有7个分量
    chosen_freq_imf = (freq >= lower) & (freq <= upper);
    imf_s1 = u2(chosen_freq_imf, i);
    Rho = corrcoef(imf_s1, P1_2_s1);
    rho(i) = Rho(1,2);
%      rho_list(i) = rho; % 将计算出的相关系数添加到向量中
end

max_corr2 = max(rho);
m = max_corr2/2;
imf_list = []; % 创建一个空矩阵用于储存所选的IMF分量
imf_list = zeros(size(u2,1),1);
for i = 1:size(u2,2)
    if rho(i) > m 
        imf_list = imf_list + u2(:,i);
    end
end
% % m = max(Rho)/2; % 设置相关系数阈值为最大相关系数的一半
% imf_list = []; % 创建一个空矩阵用于储存所选的IMF分量
% for i = 1:7 % 遍历每个IMF分量
%     if Rho(i) >= m % 如果相关系数大于等于阈值
%         imf_list = [imf_list, u2(:, i)]; % 将该IMF分量添加到矩阵中
%     end
% end




% % 计算每个 IMF 分量与泄漏信号的相关系数（s3）
% rho_list = zeros(1, 7); % 存储相关系数
% for i = 1:7
%     u2_i = u2(:, i+1);
%     u2_i_s = u2_i(chosen_freq);
%     rho = corr(u2_i_s', P1_3_s1');
%     rho_list(i) = rho;
% end

% 计算重构信号的Hankel矩阵
L = length(imf_list);
K1 = 2000;            % 选取的Hankel矩阵的列数
H = hankel(imf_list(1:K1),imf_list(K1:L));

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
imf_list_filter = U_r * S_r * V_r';

% 计算阈值
threshold = mean(H(:)); % 或者使用中值，如：median(H(:));

% 进行信号滤波
s_filter = medfilt1(imf_list_filter, 3); % 或者使用其他滤波器，如：sgolayfilt(s_recon_filter, 3, 5);
s_filter(s_filter<threshold) = 0;

%对信号进行规范化调整
s_filter_transpose = s_filter';
s_filter_col = s_filter_transpose(:,1);
% s_filter_col = s_filter_col(1:3001);

imf_list_col = imf_list(:,1);
imf_list_col = imf_list_col(1:4001);


% 计算信噪比
noise = imf_list_col- s_filter_col;

SNR = 10*log10(norm(imf_list)^2/norm(noise)^2);


% 设置时间尺度
scales = 2.^(0:4);
K=6;
%计算权重向量
N = size(imf_list,1); % 数据长度
m = length(scales); % 分段数
w = zeros(m,1);
imfs=u2;
R = 0.15;
for i=1:m % 遍历每个尺度
    seg = reshape(imf_list(1:scales(i)*floor(N/scales(i)),:)',[],scales(i)); % 将数据划分成多个子段
    med = median(seg); % 计算每个子段的中位数
    mad = median(abs(seg - med)); % 计算 MAD
    w(i) = 1/mad; % 将权重设置为1/MAD
end
w = w./sum(w); % 标准化权重

% 计算每个IMF分量在不同时间尺度下的复合多尺度熵
imf_cmse = zeros(K, length(scales));
for i = 1:K
    imf = u2(:, i);
    for j = 1:length(scales)
 %       [~, cmse] = cmse(imf, time_scales(j));
         cmse_val = cmse(imf, scales(j), 1);
        imf_cmse(i,j) = cmse_val;
    end
end

% 绘制各IMF分量在不同时间尺度下的复合多尺度熵图
for i = 1:K
    subplot(3,2,i)
    plot(log2(scales), imf_cmse(i,:), '-o', 'LineWidth', 1.2);
    xlabel('log2(时间尺度)');
    ylabel('复合多尺度熵');
    title(['第', num2str(i), '个IMF分量']);
    grid on;
end
scales = fix(log2(length(u2(:,1))/4)); % 时间尺度的个数
num_imfs = size(imf_list, 2); % IMF 分量的个数

imfs=u2;
tau = 1;
m = 6;
R = 0.15;
scales = 2.^(0:4); % 设置时间尺度

imfs_CMSE = cell(6, 1); % 用于保存每一个 IMF 分量在不同时间尺度下的 RCMSE、CMSE、MSE、MSFE、GMSE
for i = 1:6
    x = imfs(:, i); % 取出第 i 个 IMF 分量
    [~, CMSE] = Ent_MS_Plus(x, tau, m, R); % 计算其在不同时间尺度下的 CMSE 参数
    imfs_CMSE{i} = CMSE;
end
%计算平均值作为起始阈值
threshold = mean(cellfun(@(x) x(1, 1), imfs_CMSE));

%输出平均值和阈值
fprintf('Average threshold: %.4f\n', threshold);
fprintf('***********************************************************\n');
fprintf('Selection Result:\n');

%根据阈值选择 IMF 分量
selected_imfs = [];
for i = 1:6
    if imfs_CMSE{i}(1, 1) <= threshold
        fprintf('IMF-%d is selected (RCMSE: %.4f)\n', i, imfs_CMSE{i}(1, 1));
        selected_imfs = [selected_imfs i];
    else
        fprintf('IMF-%d is not selected (RCMSE: %.4f)\n', i, imfs_CMSE{i}(1, 1));
    end
end
% 计算需要选择的 IMF 分量的加权和
reconstructed_signal = sum(imfs(:, selected_imfs), 2);

% 绘制原始信号和重构信号的图像
figure;
plot(s2);
hold on;
plot(reconstructed_signal, '--');
xlabel('Time');
ylabel('Amplitude');
title('Original Signal vs. Reconstructed Signal');
legend({'Original Signal', 'Reconstructed Signal'});


