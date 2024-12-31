function [u, u_hat, omega] = VMD(s, alpha, tau, K, DC, init, tol)
% 变分模式分解
%输入和参数：
%信号-要分解的时域信号（1D）
%alpha——数据保真度约束的平衡参数
%tau-双上升的时间步长（选择0表示噪声松弛）
%K-要恢复的模式数
%DC-如果第一个模式被置于并保持为DC（0-freq），则为真
%init-0=所有omegas从0开始
%1=所有ω气开始均匀分布
%2=所有omegas随机初始化
%tol——收敛准则的容差；通常在1e-6左右
%输出：
%u-分解模式的集合
%u_hat-模式的光谱
%ω-估计模式中心频率
%K.Dragomiretskiy，D.Zosso，变分模式分解，IEEE Trans。
% 输入信号的周期和采样频率
save_T = length(s);
fs = 1/save_T;
 
% 通过镜像扩展信号
T = save_T;
f_mirror(1:T/2) = s(T/2:-1:1);
f_mirror(T/2+1:3*T/2) = s;
f_mirror(3*T/2+1:2*T) = s(T:-1:T/2+1);
f = f_mirror;
 
% 时域0到T（镜像信号的）
T = length(f);
t = (1:T)/T;
 
% Spectral Domain discretization
freqs = t-0.5-1/T;
 
% 谱域离散化
N = 500;
Alpha = alpha*ones(1,K);
f_hat = fftshift((fft(f)));
f_hat_plus = f_hat;
f_hat_plus(1:T/2) = 0;
 
% matrix keeping track of every iterant // could be discarded for mem
u_hat_plus = zeros(N, length(freqs), K);
 
% Initialization of omega_k
omega_plus = zeros(N, K);
switch init
    case 1
        for i = 1:K
            omega_plus(1,i) = (0.5/K)*(i-1);
        end
    case 2
        omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
    otherwise
        omega_plus(1,:) = 0;
end
 
% if DC mode imposed, set its omega to 0
if DC
    omega_plus(1,1) = 0;
end
 
% start with empty dual variables
lambda_hat = zeros(N, length(freqs));
 
% other inits
uDiff = tol+eps; % update step
n = 1; % loop counter
sum_uk = 0; % accumulator
%%-------迭代更新的主循环
while ( uDiff > tol &&  n < N ) % not converged and below iterations limit
     
    % update first mode accumulator
    k = 1;
    sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
     
    % update spectrum of first mode through Wiener filter of residuals
    u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
     
    % update first omega if not held at 0
    if ~DC
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
    end
     
    % update of any other mode
    for k=2:K
         
        % accumulator
        sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
         
        % mode spectrum
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
         
        % center frequencies
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
         
    end
     
    % Dual ascent
    lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
     
    % loop counter
    n = n+1;
     
    % converged yet?
    uDiff = eps;
    for i=1:K
        uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
    end
    uDiff = abs(uDiff);
     
end
 
 
%后处理和清理
 
 
% discard empty space if converged early
N = min(N,n);
omega = omega_plus(1:N,:);
 
% Signal reconstruction
u_hat = zeros(T, K);
u_hat((T/2+1):T,:) = squeeze(u_hat_plus(N,(T/2+1):T,:));
u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:)));
u_hat(1,:) = conj(u_hat(end,:));
 
u = zeros(K,length(t));
 
for k = 1:K
    u(k,:)=real(ifft(ifftshift(u_hat(:,k))));
end
 
% remove mirror part
u = u(:,T/4+1:3*T/4);
 
% recompute spectrum
clear u_hat;
for k = 1:K
    u_hat(:,k)=fftshift(fft(u(k,:)))';
end
 
end