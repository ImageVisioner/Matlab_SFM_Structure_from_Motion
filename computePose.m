% 函数：计算相机2的外参数
% 用本质矩阵的奇异值分解求解
function Rt_best = computePose(E,K,matchedpoint1,matchedpoint2)
[U,~,V] = svd(E);
W = [0 -1 0;1 0 0;0 0 1];
Z = [0 1 0;-1 0 0;0 0 0];
S = U*Z*U';
R1 = U*W*V';
R2 = U*W'*V';
t1 = U(:,3);
t2 = -U(:,3);
if det(R1) < 0
  R1 = -R1;
end

if det(R2) < 0
  R2 = -R2;
end

% 本质矩阵的奇异值分解会产生4种情况。见《计算机视觉多视图几何》第八章8.6节
Rt1 = [R1 t1];
Rt2 = [R1 t2];
Rt3 = [R2 t1];
Rt4 = [R2 t2];
Rt = {Rt1;Rt2;Rt3;Rt4};

Rt0 = [eye(3) [0 0 0]'];

% 找到正确的解
for i = 1:4
    Rt_best = findTureP(K,Rt0,Rt{i},matchedpoint1,matchedpoint2);
    if (size(Rt_best,1)~=0)
        break;
    end
end

end %end of function