function Positions=initialization(SearchAgents_no,dim,ub,lb)
% 子函数初始化种群位置, 作者邮箱：fty_mat@163.com

Positions = zeros(SearchAgents_no,dim);
for i=1:dim
    ub_i=ub(i);
    lb_i=lb(i);
    Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;%产生范围内的随机数
end
