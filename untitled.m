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