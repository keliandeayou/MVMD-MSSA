function [imf,alpha,K]=IAMSSA_MVMD(varargin)
% 子函数用于改进的麻雀搜索算法优化VMD的惩罚系数alpha和分解层数K
SearchAgents_no = cell2mat(varargin(1));   % 种群大小
Max_iteration = cell2mat(varargin(2));   % 最大迭代次数
dim = cell2mat(varargin(3));   % 变量个数
lb = cell2mat(varargin(4));   % alpha范围 K范围   下限
ub = cell2mat(varargin(5));  % 上限
try
    imf_write = cell2mat(varargin(6));
catch
    imf_write = 1;   % imf_write为0，优化结果不写入EXCEL，为1，则写入结果到EXCEL中
end

t = evalin('base', 't');   % 采样时间
s =  evalin('base', 's');    % 原始数据
fs =   evalin('base', 'fs');  % 采样频率
ST = 0.2;%预警值
PD = 0.7;%发现者的比列，剩下的是加入者
SD = 0.2;%意识到有危险麻雀的比重

PDNumber = round(SearchAgents_no*PD); %发现者数量
SDNumber = round( SearchAgents_no*SD);%意识到有危险麻雀数量
if(max(size(ub)) == 1)
    ub = ub.*ones(1,dim);
    lb = lb.*ones(1,dim);
end

%% 改进点：Cat映射初始化种群
X0=initializationNew(SearchAgents_no,dim,ub,lb);
%% 改进点：精英反向
X = initializationJY(X0,SearchAgents_no,dim,ub,lb,@(x)fun(x,s));
X = min(max(X,lb),ub);
%计算初始适应度值
fitness = zeros(1,SearchAgents_no);
for i = 1:SearchAgents_no
    fitness(i) =  fun(X(i,:),s);
end
[fitness, index]= sort(fitness);%排序
BestF = fitness(1);
WorstF = fitness(end);
GBestF = fitness(1);%全局最优适应度值
for i = 1:SearchAgents_no
    X(i,:) = X0(index(i),:);
end
Convergence_curve=zeros(Max_iteration,1);
GBestX = X(1,:);%全局最优位置
X_new = X;
gbest_array = zeros(Max_iteration,2);
for i = 1: Max_iteration
    
    %% 改进点：比例系数改进
    ST = 0.2;%预警值
    b = 0.7;%比例系数
    k = 0.2;%扰动因子
    r = b*(tan(-pi*i/(4*Max_iteration) + pi/4))-k*rand();
    %防止r为负数或者大于1的数
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
            %% 改进点：改进探索者位置更新公式
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
            %产生-1，1的随机数
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
    %边界控制
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
    %更新位置
    for j=1:SearchAgents_no
        fitness_new(j) = fun(X_new(j,:),s);
    end
    for j = 1:SearchAgents_no
        if(fitness_new(j) < GBestF)
            GBestF = fitness_new(j);
            GBestX = X_new(j,:);
        end
    end
    %% 改进点：tent扰动和柯西扰动
    Avgf = mean(fitness_new);
    x0 = rand;
    for j = 1:SearchAgents_no
        if fitness_new(j)<Avgf %柯西变异
            Temp = X_new(j,:);
            r = tan((rand() - 0.5)*pi);%柯西随机数
            Temp = Temp + r.*Temp;
            Temp(Temp>ub) = ub(Temp>ub);
            Temp(Temp<lb) = lb(Temp<lb);
            fitTemp = fun(Temp,s);
            if fitTemp<fitness_new(j)
                fitness_new(j) = fitTemp;
                X_new(j,:) = Temp;
            end
        else%Tent扰动
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
    %排序更新
    [fitness, index]= sort(fitness);%排序
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

%% 结束进化 作出适应度曲线
flag=0;   % flag为1，则添加文本水印，为0，不添加水印
figure(1)
plot(Convergence_curve,'k-.o','linewidth',1)
legend('最佳适应度')
grid on
xlabel('进化代数')
ylabel('包络熵')
title('IAMSSA进化曲线')
fig_text(flag)

figure(2)
plot(gbest_array(:,1),'ks-.','linewidth',1)
grid on
xlabel('进化代数')
ylabel('惩罚因子')
title('惩罚因子的优化过程曲线')
fig_text(flag)

figure(3)
plot(round(gbest_array(:,2)),'kv-.','linewidth',1)
grid on
xlabel('进化代数')
ylabel('分解的模态数')
title('分解模态数的优化过程曲线')
fig_text(flag)

%% 最优参数
alpha=Best_pos(1);  % 惩罚因子，也称平衡参数
K=round(Best_pos(2));  % 分解的模态数
tau = 0;            % 噪声容忍度
DC = 0;             % 无直流分量
init = 1;           % 初始化中心频率为均匀分布
tol = 1e-7;         % 收敛准则容忍度

%% 优化后的VMD分解
% u：分解模式的集合
% u_hat：模式的频谱
% omega：估计模式中心频率
[u, u_hat, omega] = MVMD(s, alpha, tau, K, DC, init, tol);
[m,~]=size(u);
m=m+1;
imf=u;
res = s'-sum(u,1);
%% 作图
% 画imf分量的时域图
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
ylabel('原始信号');
xlabel('t/s')
hold off
legend('IAMSSA_VMD分解后的时域图')

% 画imf分量的频谱图
figure(5)
for i=1:m-1
    subplot(m,1,i)
    %% FFT 变换
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
legend('IAMSSA_VMD分解后的频域图')
% 打印优化结果
fprintf_result(Best_pos,s);

% % 保存结果
if imf_write==1
    xlswrite('imf结果.xlsx',u','u','A1');
    xlswrite('imf结果.xlsx',u_hat,'uhat','A1');
    xlswrite('imf结果.xlsx',omega,'omega','A1');
end

