function X = trangulate2(K,Rt1,Rt2,matchedPoint1,matchedPoint2)
point1 = matchedPoint1.Location;
point2 = matchedPoint2.Location;
num_points = size(point1,1);
p1 = K*Rt1;
p2 = K*Rt2;
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