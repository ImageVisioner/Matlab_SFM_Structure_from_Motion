function [Rt_next,x,ind] = getRt2andX(matchPointsIndex,points,X,Rt,K)
% 得到第一组匹配（图像1和图像2）中图像2的匹配点的索引和3D坐标的对应
Index_X = cell(size(X,2),2);
for i = 1:size(X,2)
    Index_X{i,1} = matchPointsIndex{1}(i,2);
    Index_X{i,2} = X(:,i);
end

% 找到第二组匹配中图像2匹配点与第一组匹配中图像2匹配点相同的索引
flag = zeros(size(matchPointsIndex{2},1),1);
for i = 1:size(matchPointsIndex{2},1)
    for j = 1:size(matchPointsIndex{1},1)
        if (matchPointsIndex{2}(i,1) == Index_X{j,1})
            flag(i) = 1;
        end
    end
end

% 筛选出第二组匹配中已知3D点的匹配点索引，并将其与对应的三维点一一对应存储
ind_cur = matchPointsIndex{2}(flag==1,:);
ind_to_X = cell(size(ind_cur,1),2);
for i = 1:size(ind_cur,1)
    for j = 1:size(Index_X,1)
        if (Index_X{j,1}==ind_cur(i,1))
            ind_to_X{i,1} = ind_cur(i,:);
            ind_to_X{i,2} = Index_X{j,2};
        end
    end
end

% 判断DLT求解PnP的条件是否满足（条件：至少6对3D与2D对应点）
if (size(ind_to_X,1)<6)
    msgbx('点数不够，无法用DLT求解PnP!','错误!');
else
    % 构造系数矩阵 A，利用奇异值分解的方法求解[R t]
    % 公式：Ax = 0;其中：
    % A = [X(:,1)' [0 0 0 0] -u*X(:,1)';
    %      [0 0 0 0] X(:,1)' -v*X(:,1)';
    %      X(:,2)' [0 0 0 0] -u*X(:,2)';
    %      [0 0 0 0] X(:,2)' -v*X(:,2)';
    %      ...     ...     ...;         ] 总共6对点，12行
    % x = [r1 r2 r3]',而 K[R,t] = [r1';r2';r3'],r1、r2、r3是3*1的列向量
    for i = 1:6
       point_cur{i} =  points{3}(ind_to_X{i,1}(:,2),:);
    end
    A = [ind_to_X{1,2}' [0 0 0 0] -1*ind_to_X{1,2}'*point_cur{1}.Location(1,1);
         [0 0 0 0] ind_to_X{1,2}' -1*ind_to_X{1,2}'*point_cur{1}.Location(1,2);
         ind_to_X{2,2}' [0 0 0 0] -1*ind_to_X{2,2}'*point_cur{2}.Location(1,1);
         [0 0 0 0] ind_to_X{2,2}' -1*ind_to_X{2,2}'*point_cur{2}.Location(1,2);
         ind_to_X{3,2}' [0 0 0 0] -1*ind_to_X{3,2}'*point_cur{3}.Location(1,1);
         [0 0 0 0] ind_to_X{3,2}' -1*ind_to_X{3,2}'*point_cur{3}.Location(1,2);
         ind_to_X{4,2}' [0 0 0 0] -1*ind_to_X{4,2}'*point_cur{4}.Location(1,1);
         [0 0 0 0] ind_to_X{4,2}' -1*ind_to_X{4,2}'*point_cur{4}.Location(1,2);
         ind_to_X{5,2}' [0 0 0 0] -1*ind_to_X{5,2}'*point_cur{5}.Location(1,1);
         [0 0 0 0] ind_to_X{5,2}' -1*ind_to_X{5,2}'*point_cur{5}.Location(1,2);
         ind_to_X{6,2}' [0 0 0 0] -1*ind_to_X{6,2}'*point_cur{6}.Location(1,1);
         [0 0 0 0] ind_to_X{6,2}' -1*ind_to_X{6,2}'*point_cur{6}.Location(1,2);];
     [~,~,V] = svd(A);
     T_cur = V(:,end);
     T(1,:) = T_cur(1:4);
     T(2,:) = T_cur(5:8);
     T(3,:) = T_cur(9:12);
     % 利用 K是上三角矩阵的性质，对T进行QR分解
     T_to_KR = T(:,1:3);
     [R_inv,K_inv] = qr(inv(T_to_KR));
     R = inv(R_inv);
     t = K_inv*T(:,4);
     Rt_next = [R t];
     disp('Rt2 = ');
     disp(Rt_next);
     % 求解第二组匹配中未知3D点的匹配点的3D坐标
     ind_cur2 = matchPointsIndex{2}(flag==0,:);
     matchedPoint1 = points{2}(ind_cur2(:,1),:);
     matchedPoint2 = points{3}(ind_cur2(:,2),:);
     X1 = trangulate2(K,Rt,Rt_next,matchedPoint1,matchedPoint2);
     x = X1;
     ind = ind_cur2;
end

end