function Positions=initialization(SearchAgents_no,dim,ub,lb)
% �Ӻ�����ʼ����Ⱥλ��, �������䣺fty_mat@163.com

Positions = zeros(SearchAgents_no,dim);
for i=1:dim
    ub_i=ub(i);
    lb_i=lb(i);
    Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;%������Χ�ڵ������
end
