% ���������ĸ������ҵ�������������ȷ�Ľ�
function Rt_best = findTureP(K,Rt0,Rt,matchedpoint0,matchedpoint1)
% ��ȡһ��ƥ�������ά��
p0 = K*Rt0;
p1 = K*Rt;
A = [p0(3,:)*matchedpoint0(1)-p0(1,:);
     p0(3,:)*matchedpoint0(2)-p0(2,:);
     p1(3,:)*matchedpoint1(1)-p1(1,:);
     p1(3,:)*matchedpoint1(2)-p1(2,:)];
[~,~,V] = svd(A);
X = V(:,end);
X = X/X(4);

x = [X(1),X(2),X(3)]';
d0 = [0 0 1]';
d1 = Rt(3,1:3)';
O0 = [0 0 0]';
O1 = (-Rt(:,1:3)'*Rt(:,4));

% �жϹ�ʽ��
% (P-O1)'*d1>0 && (P-O2)'*d2>0
% ʽ�У�P ����ά�����꣬����Σ���������O1 �ǵ�һ��������Ĺ������꣬��������
% O2 �ǵڶ�������Ĺ������꣬��������d1��d2 �ֱ������������ķ�����������������

if ((x-O0)'*d0>0) && ((x-O1)'*d1>0)
    Rt_best = Rt;
else
    Rt_best = [];
end
end