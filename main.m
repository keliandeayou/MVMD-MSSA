%% �Ľ��Ļ�����ȸ�����㷨�Ż����ģ̬�ֽ����
    clc
    clear
    close all
    % ��ȡ����
    load('SJ.mat');
    
    fs=6000;
    % ��ȡǰ10000���ȵ��ź�
    len=6000;
    s=SJ041(1:6000,2);
    % s3=SJ031(1:6000,3);
    
    % ����ʱ��
    t = (0:len-1)/fs;0;
    
    
    %% �趨�Ľ�����ȸ�����㷨����
    popsize =25;   % ��Ⱥ��С,�ɸ���
    iter = 25;   % ����������,�ɸ���
    dim = 2;   % VMD��������
    lb = [1000 6];   % alpha��Χ K��Χ����
    ub = [3000 10];  % ����
    % %��ȸ�Ż��㷨
    % [fMin , bestX,Convergence_curve1 ] = SSA(popsize, iter,lb,ub,dim,@(x)fun(x,s));
    % %��������Ⱥ�Ż��㷨
    % [gBestScore,gBest,cg_curve]=PSO(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % %��
    % [Leader_score,Leader_pos,Convergence_curve]=WOA(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % %WSO
    % [fmin0,gbest,ccurve]=WSO(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % %MFO
    % [Best_flame_score,Best_flame_pos,Convergence_curve]=MFO(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % %
    % [Alpha_score,Alpha_pos,Convergence_curve]=GWO(popsize,iter,lb,ub,dim,@(x)fun(x,s));
    % flag=0;   % flagΪ1��������ı�ˮӡ��Ϊ0�������ˮӡ
    % fig_show(1)
    % plot(Convergence_curve1,'k-.o','linewidth',1)
    % legend('�����Ӧ��')
    % grid on
    % xlabel('��������')
    % ylabel('������')
    % title('IAMSSA��������')
    % fig_text(flag)
    % 
    % %����EXCEL
    % writetable(Convergence_curve, "1.csv")
    %% �Ľ�����ȸ�����㷨SSA�Ż�VMD����
    tic ,  % ��ʼ��ʱ
    trainlabel=1;%�Ƿ�����ѵ����1��������ѵ����0����ֱ�ӵ�ȡ�ϴ�Ѱ�ŵõ��Ľ��
    if trainlabel==1
        [~,alpha,K]=IAMSSA_MVMD(popsize,iter,dim,lb,ub,1);  % 0��ʾ������IMF��1��ʾ����IMF������
        save result alpha K 
    else
        load result alpha K 
    end
toc,  % ������ʱ
tau = 0;            % �������̶�
DC = 0;             % ��ֱ������
init = 1;           % ��ʼ������Ƶ��Ϊ���ȷֲ�
tol = 1e-7;         % ����׼�����̶�
[u, ~, ~] = VMD(s, 1109.1, tau, 8, DC, init, tol);

fs=6000;
% ��ȡǰ10000���ȵ��ź�
len=6000;
s2=SJ072(1:6000,2);
s3=SJ072(1:6000,3); 
% s3=SJ014(1:6000,3);
% ����ʱ��
t = (0:len-1)/fs;

%���׷���
s2_pre = wiener2(s2); % ���ź�s2����ȥ�봦��
s3_pre = wiener2(s3); % ���ź�s3����ȥ�봦��

fs = 6000; % ����Ƶ�ʣ���λHz
N = length(s2_pre);
S2_fft = fft(s2_pre, N)/N; % ʹ��MATLAB�е�fft�������и���Ҷ�任
S3_fft = fft(s3_pre, N)/N; % ʹ��MATLAB�е�fft�������и���Ҷ�任

f = fs*(0:(N/2))/N; % ����Ƶ����������λHz
S2_psd = 2*abs(S2_fft(1:N/2+1)).^2; % �����źŵĹ������ܶ�
S3_psd = 2*abs(S3_fft(1:N/2+1)).^2; % �����źŵĹ������ܶ�

[P_S2S3, f] = cpsd(s2_pre, s3_pre, [], [], N, fs); % �����źŵĻ����ܶ�
P_S2S3_norm = abs(P_S2S3)./sqrt(S2_psd.*S3_psd); % ��һ�������õ����׵���س̶�

omega_s1 = 1000;  % rad/s
omega_s2 = 1000;  % rad/s
omega = (2*pi*f)';
index = find((omega>=omega_s1)&(omega<=omega_s2)); % ���ҹ�������Ƶ�� ��s1����s2 ������
P_S2S3_norm_filt = P_S2S3_norm(index); % �ӹ�һ����Ļ�������ȡ��������Ƶ�� ��s1����s2 ����Ϣ



%% �����Ż����ò�������VMD�ֽ�
% load result alpha K %�������Ų���
% alpha=2000;
% K=7;
clc
clear
close all
% ��ȡ����
load('SJ.mat');

fs=6000;
% ��ȡǰ10000���ȵ��ź� 
len=6000;
s2=SJ072(1:6000,2);
s3=SJ072(1:6000,3); 

% ����ʱ��
t = (0:len-1)/fs;
% alpha=2000;
% K=7;
tau = 0;            % �������̶�
DC = 0;             % ��ֱ������
init = 1;           % ��ʼ������Ƶ��Ϊ���ȷֲ�
tol = 1e-7;         % ����׼�����̶�
data2=s2;             %dataΪ���ֽ����ݣ�������Ϊ���ֽ����ݼ��е�����һ���������������ֱ�ӱ�дһ��ѭ��ʵ���������ݵķֽ�
data3=s3;             %dataΪ���ֽ����ݣ�������Ϊ���ֽ����ݼ��е�����һ���������������ֱ�ӱ�дһ��ѭ��ʵ���������ݵķֽ�

% [u2, ~, ~] = MVMD(s2, alpha, tau, K, DC, init, tol);
% [u3, ~, ~] = MVMD(s3, alpha, tau, K, DC, init, tol);
[u2, ~, ~] = VMD(s2, 605, tau, 3, DC, init, tol);%�������Ų������зֽ�
[u3, ~, ~] = VMD(s3, 500, tau, 3, DC, init, tol);%�������Ų������зֽ�
%������� u��һ����������ÿһ�д���һ��IMF������ÿһ�д���һ��ʱ��㡣
%K1+1��IMF�ķֽ������ʵ��ѡȡ�ֽ�Ĳ���ΪK1��Ҳ����VMD�ֽ�õ���IMF��������K1��������Ҫ����������ķ�����ȥ����ѡȡǰK1��IMF���������ع���s2_vmd = sum(IMF(1:K1, :))��


%�����������ԭ�ź����ϵ��
v2 = u2';             % ת�õõ�������
v3 = u3';             % ת�õõ�������
kr2 = zeros(1, 3);   % ���ڴ洢���ϵ��������
for i = 1:3
    [Ra2, ~] = xcorr(v2(:, i), s2);
    [Rb2, lag2] = xcorr(s2);% ��������غ���Rb
    C2 = corrcoef(Ra2, Rb2);
    kr2(i) = C2(1, 2);
end
% ������״ͼ
bar(kr2);
title('��������ԭ�ź����ϵ��');
xlabel('����');
ylabel('���ϵ��');

kr3 = zeros(1, 3);   % ���ڴ洢���ϵ��������
for i = 1:3
    [Ra3, ~] = xcorr(v3(:, i), s3);
    [Rb3, lag3] = xcorr(s3);% ��������غ���Rb
    C3 = corrcoef(Ra3, Rb3);
    kr3(i) = C3(1, 2);
end
% ������״ͼ
bar(kr3);
title('��������ԭ�ź����ϵ��');
xlabel('����');
ylabel('���ϵ��');


% ѡȡ���ϵ������������ϵ��һ��ķ��������ع�
max_corr2 = max(kr2);
threshold2 = max_corr2/3;
xr2 = zeros(size(v2,1),1);

for i = 1:size(v2,2)
    if kr2(i) > threshold2
        xr2 = xr2 + v2(:,i);
    end
end

% ѡȡ���ϵ������������ϵ��һ��ķ��������ع�
max_corr3 = max(kr3);
threshold3 = max_corr3/3;
xr3 = zeros(size(v3,1),1);

for i = 1:size(v3,2)
    if kr3(i) > threshold3
        xr3 = xr3 + v3(:,i);
    end
end


% ����ع������������
disp('�ع������������');
disp(xr2');

% ����ع������������
disp('�ع������������');
disp(xr3');


%�ع��ź���ԭʼ�źŶԱ�ͼ
t2 = 1:length(s2); %����ʱ����
figure;
plot(t2,s2,'b');
hold on;
plot(t2,xr2,'r');
title('ԭʼ�źź��ع��źŶԱ�');
xlabel('ʱ��');
ylabel('��ֵ');
legend('ԭʼ�ź�','�ع��ź�');
%�����Ա�ͼ�ֿ���
t2 = 1:length(s2); %����ʱ����
figure;
subplot(2,1,1);
plot(t2,s2);
title('ԭʼ�ź�');
xlabel('ʱ��');
ylabel('��ֵ');
subplot(2,1,2);
plot(t2,xr2);
title('�ع��ź�');
xlabel('ʱ��');
ylabel('��ֵ');



%�ع��ź���ԭʼ�źŶԱ�ͼ
t3 = 1:length(s3); %����ʱ����
figure;
plot(t3,s3,'b');
hold on;
plot(t3,xr3,'r');
title('ԭʼ�źź��ع��źŶԱ�');
xlabel('ʱ��');
ylabel('��ֵ');
legend('ԭʼ�ź�','�ع��ź�');
%�����Ա�ͼ�ֿ���
t3 = 1:length(s3); %����ʱ����
figure;
subplot(2,1,1);
plot(t3,s3);
title('ԭʼ�ź�');
xlabel('ʱ��');
ylabel('��ֵ');
subplot(2,1,2);
plot(t3,xr3);
title('�ع��ź�');
xlabel('ʱ��');
ylabel('��ֵ');


x=s_filter_col;
y=s_filter_col3;
Fs=6000;
% ���������ź����ݷֱ�Ϊ s1 �� s2������Ƶ��Ϊ Fs�����ݳ��ȷֱ�Ϊ N1 �� N2
% ��Ҫ��ǰȷ�� N1 �� N2 ��Сһ��

% �����ź� s1 �� s2 �Ļ���غ���
corr = xcorr(x, y);

% ��ȡ����غ����ĳ���
M = length(corr);

% ����ʱ���ӳ٣���λ���룩
% ���� Fs = 1��ʱ���ӳٵĵ�λ���룻���� Fs ������ 1��ʱ���ӳٵĵ�λ��������
dt = (-floor(M/2):floor(M/2))/Fs;

% �ҳ�����غ��������ֵ����λ�ã���Ϊ�����ź�֮���������ֵ��ʱ���ӳ�
[max_corr, max_pos] = max(abs(corr));
delay_time = dt(max_pos);

% ������
fprintf('ʱ���ӳ�Ϊ %.2f �루���� %d ��������\n', delay_time, round(delay_time*Fs));


% %����ط���
% x=xr2;
% y=xr3;
% [coeff, lag] = xcorr(x, y);
% [~, max_index] = max(coeff);
% delay = lag(max_index);
% sampling_rate = 6000; % ���������Ϊ 6000
% delay_time = delay / sampling_rate;


% % ���ƻ����ͼ
% figure;
% plot(lag, coeff);
% title('X��Y�Ļ����');
% xlabel('ʱ��/s');
% ylabel('���ϵ��');


% �ܵ�й©�����ź��عܵ��Ĵ����ٶ�
c = 1478; % m/s
% ����ʱ���ӳټ���й©����봫������λ��
d=6.6; % m
d2 = (d+c * delay_time) / 2;







 

 
 
 %% ����uΪ������������NΪ��������Ŀ����MΪÿ�������ĳ���
%% ��ʼ�����ϵ������ÿ������֮������ϵ�����վ���Գ��ԶԳ�ȡֵ
corr_mat = corrcoef(u);
% ������ϵ������
disp('����֮������ϵ������');
disp(corr_mat);


% ��������Ծ���
figure;
imagesc(abs(corr_mat));
colorbar;
colormap(jet);
title('����֮������ϵ��');
xlabel('�������');
ylabel('�������');





kurt_vals = kurtosis(u');

kurt_thresh = 2; % ����ѡ���Ͷ�ֵ����2��ģ̬���������ع�
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
% figure('Name','Ƶ��ͼ','Color','white');
% for i = 1:K
% p=abs(fft(u(:,i))); %��fft���õ�p�����ǰ����ߵ�fft---������ 
% subplot(K,1,i);
% plot((0:len-1)*fs/len,p) %���ư�����
% xlim([0 fs/2]) %չʾ�����׵�Ƶ�Σ�����������Լ��������ѡ���Ƿ�ע��
% if i ==1
% title('Ƶ��ͼ'); xlabel('Ƶ��'); ylabel(['IMF' int2str(i)]);%int2str(i)�ǽ���ֵi���������ת����ַ���y������
% else
% xlabel('Ƶ��'); ylabel(['IMF' int2str(i)]);%int2str(i)�ǽ���ֵi���������ת����ַ���y������
% end
% end
% set(gcf,'color','w');
%% ��ͼ
% ��imf������ʱ��ͼ
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
title('VMD�ֽ���ʱ��ͼ')
% ��imf������Ƶ��ͼ
figure
for i=1:m-1
    subplot(m,1,i)
    %% FFT �任
    [cc,y_f]=plot_fft(u(i,:),fs,1);%�������ʱ���ź�u(i,:)����FFT�任
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
title('VMD�ֽ���Ƶ��ͼ')
%% δ�Ż�VMD�ֽ���
alpha=1000;  % �ͷ����ӣ�Ҳ��ƽ�����
K=5;  % �ֽ��ģ̬��
tau = 0;            % �������̶�
DC = 0;             % ��ֱ������
init = 1;           % ��ʼ������Ƶ��Ϊ���ȷֲ�
tol = 1e-7;         % ����׼�����̶�
%% δ�Ż�VMD�ֽ�
% u���ֽ�ģʽ�ļ���
% u_hat��ģʽ��Ƶ��
% omega������ģʽ����Ƶ��
[u, u_hat, omega] = VMD(s, alpha, tau, K, DC, init, tol);
%% ��ͼ
% ��imf������ʱ��ͼ
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
title('VMD�ֽ���ʱ��ͼ')
% ��imf������Ƶ��ͼ
figure
for i=1:m-1
    subplot(m,1,i)
    %% FFT �任
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
title('VMD�ֽ���Ƶ��ͼ')
tic
clc
clear all




clc
clear
close all
% ��ȡ����
load('SJ.mat');

fs=6000;
% ��ȡǰ6000���ȵ��ź�
len=6000;
s=SJ110(1:6000,2);
% ����ʱ��
t = (0:len-1)/fs;
% ����MVMD�������зֽ�
[numIMFs,IMF] = MVMD(s);
%2.����MVMD�ֽ���
vmd = matvmd(s);

%3.����MVMD�ֽ�
[Coef,v] = vmd.vmd();

%4.ѡ����ʵķ���
energy = sum(abs(Coef).^2,2);
energy_ratio = energy./sum(energy);
components = v(:,energy_ratio>0.1);

%���Ʒֽ���
figure;
subplot(2,1,1);
plot(signal);
title('Original Signal');
subplot(2,1,2);
plot(sum(Coef,1)');
title('Sum of Coefficients');

