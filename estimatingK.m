%%函数：计算内参矩阵
%通过3个消影点找无穷远平面，然后自标定求相机内参数

function [K] = estimatingK(img_rgb)

% 计算消影点
vps = detectVP(img_rgb);


% 绝对二次曲线图像是二次曲线：w = (KK^T)^-1
% W = [w1,  0, w2;
%       0, w1, w3;
%      w2, w3, w4]
% w = [w1,w2,w3,w4]'
% Aw = 0
% 公式参考《计算机视觉中的多视图几何》第七章的7.7节
A = [vps(1, 1)*vps(2,1)+vps(1,2)*vps(2,2), vps(1,1)+vps(2,1), vps(1,2)+vps(2,2), 1;
     vps(2, 1)*vps(3,1)+vps(2,2)*vps(3,2), vps(2,1)+vps(3,1), vps(2,2)+vps(3,2), 1;
     vps(3, 1)*vps(1,1)+vps(3,2)*vps(1,2), vps(3,1)+vps(1,1), vps(3,2)+vps(1,2), 1];
[~, ~, v] = svd(A);
w = v(:, end);
W = [w(1), 0, w(2); 0, w(1), w(3); w(2), w(3), w(4)];
K = inv(chol(W));
K = K / K(3, 3);

end