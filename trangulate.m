% 函数：三角化求匹配点的三维坐标
function X = trangulate(K,Rt,matchedPoint1,matchedPoint2)
point1 = matchedPoint1.Location;
point2 = matchedPoint2.Location;
num_points = size(point1,1);
p1 = K*[eye(3) [0 0 0]'];
p2 = K*Rt;
for i = 1:num_points
    A = [p1(3,:)*point1(i,1)-p1(1,:);
        p1(3,:)*point1(i,2)-p1(2,:);
        p2(3,:)*point2(i,1)-p2(1,:);
        p2(3,:)*point2(i,2)-p2(2,:)];
    [~,~,V] = svd(A);
    X(:,i) = V(:,end);
    X(:,i) = X(:,i)/X(end,i);
end
end