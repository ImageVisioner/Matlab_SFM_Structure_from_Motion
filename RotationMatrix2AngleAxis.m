function angle_axis = RotationMatrix2AngleAxis(R)

% 此函数是将旋转矩阵转变为旋转向量的形式，则旋转矩阵R将变为旋转向量r
% r = (r1,r2,r3),theta是向量r的模长

% 记旋转向量r的单位向量r_unit = (rx,ry,rz)
% sin(theta)*[0 -rz ry; rz 0 -rx; -ry rx 0] = (R - R')/2
% angle_axis(1) = 2*sin(theta)*rx = R(3,2)-R'(2,3);
% angle_axis(2) = 2*sin(theta)*ry = R(1,3)-R'(3,1);
% angle_axis(3) = 2*sin(theta)*rz = R(2,1)-R'(1,2);                    
angle_axis(1) = R(1+2, 1+1) - R(1+1, 1+2);
angle_axis(2) = R(1+0, 1+2) - R(1+2, 1+0);
angle_axis(3) = R(1+1, 1+0) - R(1+0, 1+1);


% cos(theta) = (trace(R)-1)/2; 注：trace为matlab自带求迹的函数
% cos(theta)应在[-1,1]之间，故对取值进行限制
% trace(R) = R(1,1)+R(2,2)+R(3,3)
costheta = min(max((R(1+0, 1+0) + R(1+1, 1+1) + R(1+2, 1+2) - 1.0) / 2.0,  -1.0), 1.0);

% 4*sin(theta)^2*(rx^2+ry^2+rz^2) = angle_axis(1)^2 + angle_axis(2)^2 + angle_axis(2)^2
% 因为r_unit为单位向量，故：rx^2+ry^2+rz^2 = 1;
% sin(theta) = sqrt(angle_axis(1)^2 + angle_axis(2)^2 + angle_axis(2)^2)/2
% sin(theta)应在[-1,1]之间，sqrt()函数所求值不会小于0，故限制上限
sintheta = min(sqrt(angle_axis(1) * angle_axis(1) + angle_axis(2) * angle_axis(2) + angle_axis(3) * angle_axis(3)) / 2.0, 1.0);

% 使用反正切求theta
theta = atan2(sintheta, costheta);

% r1 = rx*theta = angle_axis(1)*theta/(2*sin(theta))
% r2 = ry*theta = angle_axis(2)*theta/(2*sin(theta))
% r3 = rz*theta = angle_axis(3)*theta/(2*sin(theta))
% 第一种情况：sin(theta)足够大，上式的计算不需要改变
% 因为若是sin(theta)太小,上式求得的值就会太大，导致数据溢出
% 所以，需要对sin(theta)的值做一个阈值检验
kThreshold = 1e-12;
if ((sintheta > kThreshold) || (sintheta < -kThreshold))
    r = theta / (2.0 * sintheta);
    angle_axis = angle_axis * r;
    return;
end

% 第二种情况，theta接近于0（表明theta不足够大），此时sin(theta)近似等于theta
% 所以可以用theta的值代替sin(theta)
if (costheta > 0.0)
    angle_axis = angle_axis * 0.5;
    return;
end

% 第三种情况，theta接近于pi,此时sin(theta)接近于0，theta不近似等于sin(theta)
% 此时情况较为复杂，运用下面的计算解决
inv_one_minus_costheta = 1.0 / (1.0 - costheta);
for i=1:3
    angle_axis(i) = theta * sqrt((R(i, i) - costheta) * inv_one_minus_costheta);
    if (((sintheta < 0.0) && (angle_axis(i) > 0.0)) || ((sintheta > 0.0) && (angle_axis(i) < 0.0)))
        angle_axis(i) = -angle_axis(i);
    end
end
