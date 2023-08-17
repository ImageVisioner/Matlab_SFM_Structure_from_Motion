% 函数：从四个解中找到满足条件的正确的解
function Rt_best = findTureP(K,Rt0,Rt,matchedpoint0,matchedpoint1)
% 任取一对匹配点求三维点
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

% 判断公式：
% (P-O1)'*d1>0 && (P-O2)'*d2>0
% 式中：P 是三维点坐标，非齐次，列向量；O1 是第一个摄像机的光心坐标，列向量；
% O2 是第二个相机的光心坐标，列向量；d1、d2 分别是两相机朝向的方向向量，列向量；

if ((x-O0)'*d0>0) && ((x-O1)'*d1>0)
    Rt_best = Rt;
else
    Rt_best = [];
end
end