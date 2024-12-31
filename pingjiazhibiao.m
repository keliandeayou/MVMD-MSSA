N = 6000;  % 分量的长度
fs =6000;    % 采样率，这里假设为1
%%将IMF从时域转到频域
% 假设 u1 到 u7 是 6000x1 的数值矩阵

umatrix = [u1, u2, u3, u4, u5, u6, u7]; % 将 u1 到 u7 放入一个矩阵

for i = 1:size(umatrix, 2)
    u = umatrix(:, i);
    n = length(u);
    num_points = 2^nextpow2(n);
    u_padded = [u; zeros(num_points-n, 1)];
    u_frequency = fft(u_padded);
    u_spectrum = abs(u_frequency);
    
    % 可选：绘制频谱图
    sampling_rate = 1; % 根据实际情况替换采样率
    frequency_axis = (0:(num_points/2))*(sampling_rate/num_points);
    
    subplot(size(umatrix, 2), 1, i);
    plot(frequency_axis, u_spectrum(1:num_points/2+1));
    xlabel('频率（Hz）');
    ylabel('幅度');
    title(sprintf('u%d 的频谱图', i));
end
%%
N = 6000;  % 分量的长度
fs =6000;    % 采样率，这里假设为1
spectrum_u1 = abs(fftshift(fft(u1)));
spectrum_u2 = abs(fftshift(fft(u2)));
spectrum_u3 = abs(fftshift(fft(u3)));
spectrum_u4 = abs(fftshift(fft(u4)));
spectrum_u5 = abs(fftshift(fft(u5)));
spectrum_u6 = abs(fftshift(fft(u6)));
spectrum_u7 = abs(fftshift(fft(u7)));
spectrum_u8 = abs(fftshift(fft(u8)));
A1 = spectrum_u1;
A2 = spectrum_u2;
A3 = spectrum_u3;
A4 = spectrum_u4;
A5 = spectrum_u5;
A6 = spectrum_u6;
A7 = spectrum_u7;
A8 = spectrum_u8;
%%计算中心频率
% U1 = fft(u1);
% U2 = fft(u2);
% U3 = fft(u3);
% U4 = fft(u4);
% U5 = fft(u5);
% U6 = fft(u6);
% A1 = abs(U1);
% A2 = abs(U2);
% A3 = abs(U3);
% A4 = abs(U4);
% A5 = abs(U5);
% A6 = abs(U6);

[maxA1, idx1] = max(A1);
[maxA2, idx2] = max(A2);
[maxA3, idx3] = max(A3);
[maxA4, idx4] = max(A4);
[maxA5, idx5] = max(A5);
[maxA6, idx6] = max(A6);
[maxA7, idx7] = max(A7);
[maxA8, idx8] = max(A8);

f_resolution = fs / N;
f1 = (idx1 - 1) * f_resolution;
f2 = (idx2 - 1) * f_resolution;
f3 = (idx3 - 1) * f_resolution;
f4 = (idx4 - 1) * f_resolution;
f5 = (idx5 - 1) * f_resolution;
f6 = (idx6 - 1) * f_resolution;
f7 = (idx7 - 1) * f_resolution;
f8 = (idx8 - 1) * f_resolution;

%%计算带宽

half_energy1 = 0.5 * sum(A1) / 2;   % 一半的能量
half_energy2 = 0.5 * sum(A2) / 2;
half_energy3 = 0.5 * sum(A3) / 2;
half_energy4 = 0.5 * sum(A4) / 2;
half_energy5 = 0.5 * sum(A5) / 2;
half_energy6 = 0.5 * sum(A6) / 2;
half_energy7 = 0.5 * sum(A7) / 2;
half_energy8 = 0.5 * sum(A7) / 2;
% 找到频谱中超过一半能量的最大索引位置
idx1 = find(cumsum(A1) >= half_energy1, 1);
idx2 = find(cumsum(A2) >= half_energy2, 1);
idx3 = find(cumsum(A3) >= half_energy3, 1);
idx4 = find(cumsum(A4) >= half_energy4, 1);
idx5 = find(cumsum(A5) >= half_energy5, 1);
idx6 = find(cumsum(A6) >= half_energy6, 1);
idx7 = find(cumsum(A7) >= half_energy7, 1);
idx8 = find(cumsum(A8) >= half_energy8, 1);
% 根据采样率和频率分辨率计算带宽
f_resolution = fs / N;

bandwidth1 = (idx1 - 1) * f_resolution;
bandwidth2 = (idx2 - 1) * f_resolution;
bandwidth3 = (idx3 - 1) * f_resolution;
bandwidth4 = (idx4 - 1) * f_resolution;
bandwidth5 = (idx5 - 1) * f_resolution;
bandwidth6 = (idx6 - 1) * f_resolution;
bandwidth7 = (idx7 - 1) * f_resolution;
bandwidth8 = (idx8 - 1) * f_resolution;
%%计算频率范围

lower_freq1 = f1 - bandwidth1 / 2;
upper_freq1 = f1 + bandwidth1 / 2;

lower_freq2 = f2 - bandwidth2 / 2;
upper_freq2 = f2 + bandwidth2 / 2;

lower_freq3 = f3 - bandwidth3 / 2;
upper_freq3 = f3 + bandwidth3 / 2;

lower_freq4 = f4 - bandwidth4 / 2;
upper_freq4 = f4 + bandwidth4 / 2;

lower_freq5 = f5 - bandwidth5 / 2;
upper_freq5 = f5 + bandwidth5 / 2;

lower_freq6 = f6 - bandwidth6 / 2;
upper_freq6 = f6 + bandwidth6 / 2;

lower_freq7 = f7 - bandwidth7 / 2;
upper_freq7 = f7 + bandwidth7 / 2;
% 
 lower_freq8 = f8 - bandwidth8 / 2;
 upper_freq8 = f8 + bandwidth8 / 2;

%%六个分量频率范围
%第一列数据(K=6)
%分量1    1552.5-4447.5
%分量2    1212.5-3613.5
%分量3    869.5-2614.5
%分量4    602 -1812
%分量5    570.5-1653.5
%分量6    346-994
%分量7    346-994
%第二列数据(K=6)
%分量1    1556-4444
%分量2    1244-3710
%分量3    1014.5-2861.5
%分量4    610 -1886
%分量5    617.5-1774.5
%分量6    185.5-646.5
%分量7    185.5-646.5

%第一列数据(K=7)
%分量1    1550-4450
%分量2    1198.5-3627.5
%分量3    1035.5-3038.5
%分量4    736 -2222
%分量5    577-1755
%分量6    533.5-1596.5
%分量7    351-989
%第二列数据(K=7)
%分量1    1554-4446
%分量2    1238.5-3715.5
%分量3    1048-3074
%分量4    753.5 -2212.5
%分量5    624.5-1871.5
%分量6    560-1664
%分量7    161.5-478.5


%第一列数据(K=8)SJ031
%分量1    1507-4393
%分量2    1208.5-3641.5
%分量3    1065.5-3122.5
%分量4    823 -2531
%分量5    646.5-1993.5
%分量6    496-1516
%分量7    304.5-949.5
%分量8    173.5-438.5
%第二列数据(K=8)SJ031
%分量1    1495.5-4404.5
%分量2    1273.5-3858.5
%分量3    1138-3310
%分量4    899 -2699
%分量5    743.5-2220.5
%分量6    525.5-1616.5
%分量7    373-1069
%分量8    64-294
%%MSSA中的P就是Hankel矩阵的秩数(0-6000)













%分量1    
%分量2    
%分量3    
%分量4    
%分量5    
%分量6    
%分量7


%%观察MSSA的去噪效果
% 假设原始信号为 u1，经过 MSSA 滤波后的信号为 DATA_f
%IMF1过滤前后比较
% 观察时域波形
figure;
subplot(2,1,1);
plot(u1);
title('原始信号');
subplot(2,1,2);
plot(DATA_f1);
title('滤波后的信号');

% 分析频谱
Fs = 6000; % 假设采样率为 1000 Hz
N = length(spectrum_u1);
f = (-Fs/2) : (Fs/N) : (Fs/2 - Fs/N);

U1 = spectrum_u1;
Spectrum_DATA_f1 = abs(fftshift(fft(DATA_f1)));

figure;
subplot(2,1,1);
plot(f, U1);
title('原始信号频谱');
subplot(2,1,2);
plot(f, Spectrum_DATA_f1);
title('滤波后的信号频谱');

% 计算噪声功率
noise_power_spectrum_u1 = mean(spectrum_u1.^2);
noise_power_DATA_f = mean(DATA_f1.^2);

fprintf('原始信号的噪声功率：%f\n', noise_power_spectrum_u1);
fprintf('滤波后的信号的噪声功率：%f\n', noise_power_DATA_f);



%%求分量滤波前后的差分数据图
% 假设有原始数据 original_data 和滤波后数据 filtered_data
original_data=u1;
filtered_data=DATA_f1;
% 计算差分数据1
diff_data = original_data - filtered_data;

% 绘制差分数据图
figure;
plot(diff_data);
xlabel('Time');
ylabel('Difference');
title('Difference Plot: Before and After Filtering');
grid on;

% 计算滤波前后信号的频谱
spectrum_original = abs(fft(original_data));
spectrum_filtered = abs(fft(filtered_data));

% 绘制频谱图
figure;
subplot(2,1,1);
plot(spectrum_original);
xlabel('频率');
ylabel('幅度');
title('滤波前的信号频谱');

subplot(2,1,2);
plot(spectrum_filtered);
xlabel('频率');
ylabel('幅度');
title('滤波后的信号频谱');





if meth==1;    
%[U,S,V] = svds(M,1); Mout = U(:,1:P)*S(1:P,1:P)*V(:,1:P)'; end;
  [U,S,V] = svds(M,1); Mout = U*S*V; end;
H = hankel(DATA(1:2000),DATA(2000:6000));

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
Mout = U_r * S_r * V_r';

if meth==2;                           Mout = rqrd(M,P); end


%%互相关分析
x=s22;
y=s33;
x=reconstructed_signal22;
y=reconstructed_signal33;

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
dt =dt';
% 绘制互相关图
figure;
plot(dt, corr);
xlabel('延迟（秒）');
ylabel('互相关函数');
title('两个重构信号的互相关图');
%dt =dt'; 
% 管道泄漏的声信号沿管道的传播速度
c = 1478; % m/s
% 根据时间延迟计算泄漏点距离传感器的位置
d=6.6; % m
d2 = (d+c * delay_time) / 2;

%%计算信噪比
% 假设你有信号数据 s2 和噪声数据 nn
s2=reconstructed_signal2;
% 计算信号功率 Ps
Ps = (1 / length(s2)) * sum(s2.^2);
% 计算噪声功率 Pn
Pn = (1 / length(n)) * sum(n.^2);
% 计算信噪比 SNR
SNR = 10 * log10(Ps / Pn);
% 显示信噪比结果
disp(['信噪比 SNR = ', num2str(SNR), ' dB']);



% 计算原始信号和重构信号的平方和
signal_power = sum(s3.^2);
error_power = sum((s3 - reconstructed_signal2222).^2);
% 计算SNR0
SNR0 = 10 * log10(signal_power / error_power);

% 2:1.9046
% 22:16.4117
% 222：1.4508；
% 2222：5.9790；
% 进行互相关分析
xc = xcorr(s2, s3);

% 找到互相关函数的峰值和位置
[max_peak, max_index] = max(xc);

% 计算峰值相对于零延迟的偏移量
zero_lag = length(s2);
peak_offset = max_index - zero_lag;

% 显示相关峰值
fprintf('相关峰值: %.2f\n', max_peak);
fprintf('相关峰值偏移量: %d\n', peak_offset);
