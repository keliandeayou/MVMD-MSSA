
%% 改进的基于麻雀搜索算法优化变分模态分解参数
clc
clear
close all
% 读取数据
load('SJ80.mat');

fs=6000;
% 读取前10000长度的信号
len=10000;
s=SJ80(1:10000,3);
%s=X118_DE_time(1:len);
% 采样时间
t = (0:len-1)/fs;


imf = cell(2, 5);  % 用于存储 VMD 分解后的 IMFs
for i = 1:2 % 对两路传感器进行分解
    [imf{i,:}, ~] = VMD(s(:,i), 5, 50, 1e-6);
end

r = cell(2, K+1); % 用于存储相关系数

for i = 1:2 % 对两路传感器进行计算
    % 计算原始信号与低频分量的相关系数
    r{i,1} = corrcoef(s(:,i), imf{i,K+1});
    r{i,1} = abs(r{i,1}(1,2));
    % 计算其他分量（IMF）与原始信号的相关系数
    for j = 1:K
        r{i,j+1} = corrcoef(s(:,i), imf{i,j});
        r{i,j+1} = abs(r{i,j+1}(1,2));
    end
end

maxR = max(cell2mat(r), [], 'all'); % 找到所有相关系数中的最大值
thresh = maxR / 2; % 设定阈值为最大值的一半
recons = cell(2, K); % 用于存储选中的分量

for i = 1:2 % 对两路传感器进行重构
    % 找到相关系数大于阈值的分量
    idx = find( cell2mat(r(i,2:end)) > thresh ) + 1;
    % 对选中的分量进行重构
    recons{i,:} = zeros(size(s,1), length(idx));
    for j = 1:length(idx)
        recons{i,j} = imf{i,idx(j)};
    end
    recons{i,end+1} = imf{i,K+1}; % 将低频分量也加入重构中
    recons{i} = sum(cat(3, recons{i,:}), 3);
end