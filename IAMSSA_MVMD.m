function [imf,alpha,K]=IAMSSA_MVMD(varargin)
% �Ӻ������ڸĽ�����ȸ�����㷨�Ż�VMD�ĳͷ�ϵ��alpha�ͷֽ����K
SearchAgents_no = cell2mat(varargin(1));   % ��Ⱥ��С
Max_iteration = cell2mat(varargin(2));   % ����������
dim = cell2mat(varargin(3));   % ��������
lb = cell2mat(varargin(4));   % alpha��Χ K��Χ   ����
ub = cell2mat(varargin(5));  % ����
try
    imf_write = cell2mat(varargin(6));
catch
    imf_write = 1;   % imf_writeΪ0���Ż������д��EXCEL��Ϊ1����д������EXCEL��
end

t = evalin('base', 't');   % ����ʱ��
s =  evalin('base', 's');    % ԭʼ����
fs =   evalin('base', 'fs');  % ����Ƶ��
ST = 0.2;%Ԥ��ֵ
PD = 0.7;%�����ߵı��У�ʣ�µ��Ǽ�����
SD = 0.2;%��ʶ����Σ����ȸ�ı���

PDNumber = round(SearchAgents_no*PD); %����������
SDNumber = round( SearchAgents_no*SD);%��ʶ����Σ����ȸ����
if(max(size(ub)) == 1)
    ub = ub.*ones(1,dim);
    lb = lb.*ones(1,dim);
end

%% �Ľ��㣺Catӳ���ʼ����Ⱥ
X0=initializationNew(SearchAgents_no,dim,ub,lb);
%% �Ľ��㣺��Ӣ����
X = initializationJY(X0,SearchAgents_no,dim,ub,lb,@(x)fun(x,s));
X = min(max(X,lb),ub);
%�����ʼ��Ӧ��ֵ
fitness = zeros(1,SearchAgents_no);
for i = 1:SearchAgents_no
    fitness(i) =  fun(X(i,:),s);
end
[fitness, index]= sort(fitness);%����
BestF = fitness(1);
WorstF = fitness(end);
GBestF = fitness(1);%ȫ��������Ӧ��ֵ
for i = 1:SearchAgents_no
    X(i,:) = X0(index(i),:);
end
Convergence_curve=zeros(Max_iteration,1);
GBestX = X(1,:);%ȫ������λ��
X_new = X;
gbest_array = zeros(Max_iteration,2);
for i = 1: Max_iteration
    
    %% �Ľ��㣺����ϵ���Ľ�
    ST = 0.2;%Ԥ��ֵ
    b = 0.7;%����ϵ��
    k = 0.2;%�Ŷ�����
    r = b*(tan(-pi*i/(4*Max_iteration) + pi/4))-k*rand();
    %��ֹrΪ�������ߴ���1����
    if r<=0
        r=0.1;
    end
    if r>1
        r=0.1;
    end
    PDNumber = round(SearchAgents_no*r);
    SDNumber = round(SearchAgents_no*(1-r));
    
    BestF = fitness(1);
    WorstF = fitness(end);
    R2 = rand(1);
    for j = 1:PDNumber
        if(R2<ST)
            %% �Ľ��㣺�Ľ�̽����λ�ø��¹�ʽ
            X_new(j,:) = X(j,:).*(2/exp(4*i./(rand()*Max_iteration)^2));
        else
            X_new(j,:) = X(j,:) + randn()*ones(1,dim);
        end
    end
    for j = PDNumber+1:SearchAgents_no
        %        if(j>(SearchAgents_no/2))
        if(j>(SearchAgents_no - PDNumber)/2 + PDNumber)
            X_new(j,:)= randn().*exp((X(end,:) - X(j,:))/j^2);
        else
            %����-1��1�������
            A = ones(1,dim);
            for a = 1:dim
                if(rand()>0.5)
                    A(a) = -1;
                end
            end
            AA = A'*inv(A*A');
            X_new(j,:)= X(1,:) + abs(X(j,:) - X(1,:)).*AA';
        end
    end
    Temp = randperm(SearchAgents_no);
    SDchooseIndex = Temp(1:SDNumber);
    for j = 1:SDNumber
        if(fitness(SDchooseIndex(j))>BestF)
            X_new(SDchooseIndex(j),:) = X(1,:) + randn().*abs(X(SDchooseIndex(j),:) - X(1,:));
        elseif(fitness(SDchooseIndex(j))== BestF)
            K = 2*rand() -1;
            X_new(SDchooseIndex(j),:) = X(SDchooseIndex(j),:) + K.*(abs( X(SDchooseIndex(j),:) - X(end,:))./(fitness(SDchooseIndex(j)) - fitness(end) + 10^-8));
        end
    end
    %�߽����
    for j = 1:SearchAgents_no
        for a = 1: dim
            if(X_new(j,a)>ub(a))
                X_new(j,a) =ub(a);
            end
            if(X_new(j,a)<lb(a))
                X_new(j,a) =lb(a);
            end
        end
    end
    %����λ��
    for j=1:SearchAgents_no
        fitness_new(j) = fun(X_new(j,:),s);
    end
    for j = 1:SearchAgents_no
        if(fitness_new(j) < GBestF)
            GBestF = fitness_new(j);
            GBestX = X_new(j,:);
        end
    end
    %% �Ľ��㣺tent�Ŷ��Ϳ����Ŷ�
    Avgf = mean(fitness_new);
    x0 = rand;
    for j = 1:SearchAgents_no
        if fitness_new(j)<Avgf %��������
            Temp = X_new(j,:);
            r = tan((rand() - 0.5)*pi);%���������
            Temp = Temp + r.*Temp;
            Temp(Temp>ub) = ub(Temp>ub);
            Temp(Temp<lb) = lb(Temp<lb);
            fitTemp = fun(Temp,s);
            if fitTemp<fitness_new(j)
                fitness_new(j) = fitTemp;
                X_new(j,:) = Temp;
            end
        else%Tent�Ŷ�
            if x0<0.5
                tentV = 2*x0+rand/SearchAgents_no;
            else
                tentV = 2*(1-x0)+rand/SearchAgents_no;
            end
            Temp = X_new(j,:);
            Temp = Temp + tentV.*Temp;
            Temp(Temp>ub) = ub(Temp>ub);
            Temp(Temp<lb) = lb(Temp<lb);
            fitTemp = fun(Temp,s);
            if fitTemp<fitness_new(j)
                fitness_new(j) = fitTemp;
                X_new(j,:) = Temp;
            end
            x0 = tentV;
        end
    end
    
    X = X_new;
    fitness = fitness_new;
    %�������
    [fitness, index]= sort(fitness);%����
    BestF = fitness(1);
    WorstF = fitness(end);
    for j = 1:SearchAgents_no
        X(j,:) = X(index(j),:);
    end
    disp(['current iteration is: ',num2str(i), ', best fitness is: ', num2str(GBestF)])
    gbest_array(i,:) = GBestX;
    Convergence_curve(i) = GBestF;
end
Best_pos =GBestX;

%% �������� ������Ӧ������
flag=0;   % flagΪ1��������ı�ˮӡ��Ϊ0�������ˮӡ
figure(1)
plot(Convergence_curve,'k-.o','linewidth',1)
legend('�����Ӧ��')
grid on
xlabel('��������')
ylabel('������')
title('IAMSSA��������')
fig_text(flag)

figure(2)
plot(gbest_array(:,1),'ks-.','linewidth',1)
grid on
xlabel('��������')
ylabel('�ͷ�����')
title('�ͷ����ӵ��Ż���������')
fig_text(flag)

figure(3)
plot(round(gbest_array(:,2)),'kv-.','linewidth',1)
grid on
xlabel('��������')
ylabel('�ֽ��ģ̬��')
title('�ֽ�ģ̬�����Ż���������')
fig_text(flag)

%% ���Ų���
alpha=Best_pos(1);  % �ͷ����ӣ�Ҳ��ƽ�����
K=round(Best_pos(2));  % �ֽ��ģ̬��
tau = 0;            % �������̶�
DC = 0;             % ��ֱ������
init = 1;           % ��ʼ������Ƶ��Ϊ���ȷֲ�
tol = 1e-7;         % ����׼�����̶�

%% �Ż����VMD�ֽ�
% u���ֽ�ģʽ�ļ���
% u_hat��ģʽ��Ƶ��
% omega������ģʽ����Ƶ��
[u, u_hat, omega] = MVMD(s, alpha, tau, K, DC, init, tol);
[m,~]=size(u);
m=m+1;
imf=u;
res = s'-sum(u,1);
%% ��ͼ
% ��imf������ʱ��ͼ
figure(4)
for i=1:K
    subplot(K+2,1,i)
    plot(t,imf(i,:),'b-','linewidth',1);hold on
    ylabel(['IMF',num2str(i)]);
end
res = s'-sum(imf,1);hold on
subplot(K+2,1,K+1)
plot(t,res,'b-','linewidth',1)
ylabel('Res');
xlabel('t/s');
hold on
subplot(K+2,1,K+2)
plot(t,s','b-','linewidth',1)
ylabel('ԭʼ�ź�');
xlabel('t/s')
hold off
legend('IAMSSA_VMD�ֽ���ʱ��ͼ')

% ��imf������Ƶ��ͼ
figure(5)
for i=1:m-1
    subplot(m,1,i)
    %% FFT �任
    [cc,y_f]=plot_fft(u(i,:),fs,1);
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
legend('IAMSSA_VMD�ֽ���Ƶ��ͼ')
% ��ӡ�Ż����
fprintf_result(Best_pos,s);

% % ������
if imf_write==1
    xlswrite('imf���.xlsx',u','u','A1');
    xlswrite('imf���.xlsx',u_hat,'uhat','A1');
    xlswrite('imf���.xlsx',omega,'omega','A1');
end

