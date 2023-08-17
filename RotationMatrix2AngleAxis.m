function angle_axis = RotationMatrix2AngleAxis(R)

% �˺����ǽ���ת����ת��Ϊ��ת��������ʽ������ת����R����Ϊ��ת����r
% r = (r1,r2,r3),theta������r��ģ��

% ����ת����r�ĵ�λ����r_unit = (rx,ry,rz)
% sin(theta)*[0 -rz ry; rz 0 -rx; -ry rx 0] = (R - R')/2
% angle_axis(1) = 2*sin(theta)*rx = R(3,2)-R'(2,3);
% angle_axis(2) = 2*sin(theta)*ry = R(1,3)-R'(3,1);
% angle_axis(3) = 2*sin(theta)*rz = R(2,1)-R'(1,2);                    
angle_axis(1) = R(1+2, 1+1) - R(1+1, 1+2);
angle_axis(2) = R(1+0, 1+2) - R(1+2, 1+0);
angle_axis(3) = R(1+1, 1+0) - R(1+0, 1+1);


% cos(theta) = (trace(R)-1)/2; ע��traceΪmatlab�Դ��󼣵ĺ���
% cos(theta)Ӧ��[-1,1]֮�䣬�ʶ�ȡֵ��������
% trace(R) = R(1,1)+R(2,2)+R(3,3)
costheta = min(max((R(1+0, 1+0) + R(1+1, 1+1) + R(1+2, 1+2) - 1.0) / 2.0,  -1.0), 1.0);

% 4*sin(theta)^2*(rx^2+ry^2+rz^2) = angle_axis(1)^2 + angle_axis(2)^2 + angle_axis(2)^2
% ��Ϊr_unitΪ��λ�������ʣ�rx^2+ry^2+rz^2 = 1;
% sin(theta) = sqrt(angle_axis(1)^2 + angle_axis(2)^2 + angle_axis(2)^2)/2
% sin(theta)Ӧ��[-1,1]֮�䣬sqrt()��������ֵ����С��0������������
sintheta = min(sqrt(angle_axis(1) * angle_axis(1) + angle_axis(2) * angle_axis(2) + angle_axis(3) * angle_axis(3)) / 2.0, 1.0);

% ʹ�÷�������theta
theta = atan2(sintheta, costheta);

% r1 = rx*theta = angle_axis(1)*theta/(2*sin(theta))
% r2 = ry*theta = angle_axis(2)*theta/(2*sin(theta))
% r3 = rz*theta = angle_axis(3)*theta/(2*sin(theta))
% ��һ�������sin(theta)�㹻����ʽ�ļ��㲻��Ҫ�ı�
% ��Ϊ����sin(theta)̫С,��ʽ��õ�ֵ�ͻ�̫�󣬵����������
% ���ԣ���Ҫ��sin(theta)��ֵ��һ����ֵ����
kThreshold = 1e-12;
if ((sintheta > kThreshold) || (sintheta < -kThreshold))
    r = theta / (2.0 * sintheta);
    angle_axis = angle_axis * r;
    return;
end

% �ڶ��������theta�ӽ���0������theta���㹻�󣩣���ʱsin(theta)���Ƶ���theta
% ���Կ�����theta��ֵ����sin(theta)
if (costheta > 0.0)
    angle_axis = angle_axis * 0.5;
    return;
end

% �����������theta�ӽ���pi,��ʱsin(theta)�ӽ���0��theta�����Ƶ���sin(theta)
% ��ʱ�����Ϊ���ӣ���������ļ�����
inv_one_minus_costheta = 1.0 / (1.0 - costheta);
for i=1:3
    angle_axis(i) = theta * sqrt((R(i, i) - costheta) * inv_one_minus_costheta);
    if (((sintheta < 0.0) && (angle_axis(i) > 0.0)) || ((sintheta > 0.0) && (angle_axis(i) < 0.0)))
        angle_axis(i) = -angle_axis(i);
    end
end
