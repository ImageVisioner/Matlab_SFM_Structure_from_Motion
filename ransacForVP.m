%%函数：估计消失点
%使用RANSAC拟合直线
%输入：
%   lines:([x1 x2 y1 y2 theta r])
%   clusters:(1,3)细胞阵列
%   cluster_id:1 ~ 3
%   flag:标志，真或假
%输出：
% vp:(2,1)消失点的x和y
% max_inliers: int，最大inliers数

function [vp,lines_inlier_best, lines_outlier_best] = ransacForVP(lines,clusters, cluster_id)
%RANSAC 设置
eps = 0.1;
iter_num = 30000;
min_num = 10;

%初始的返回值
line_ids = clusters{cluster_id};
num = size(line_ids,2);
max_inliers = 0;
vp = zeros(3,1);
lines_inlier_best = [];
lines_outlier_best = [];

%集群中必须有足够的直线
if num < min_num
    return;
end
for i = 1:iter_num
    lines_inlier_cur = [];
    lines_outlier_cur = [];
    cur_inliers = 0;
    
    %每次迭代随机采样2个线段
    samples = randsample(num, 2);
    id1 = line_ids(samples(1));
    id2 = line_ids(samples(2));
    
    %射影坐标系下：过两点的直线等于两点齐次向量的叉乘：I=(a1)X(a2) 
    %线段1
    pt1_1 = [lines(id1, 1), lines(id1, 3), 1];
    pt1_2 = [lines(id1, 2), lines(id1, 4), 1];
    line1 = cross(pt1_1, pt1_2);
    
    %线段2
    pt2_1 = [lines(id2, 1), lines(id2, 3), 1];
    pt2_2 = [lines(id2, 2), lines(id2, 4), 1];
    line2 = cross(pt2_1, pt2_2);
    
    %射影坐标系下：两直线的交点等于两直线齐次向量的叉乘：x=(I1)X(I2)
    %计算消失点
    vp_can = cross(line1, line2);
    
    %将齐次点最后一个元素变为1，则前两个元素为欧式空间的坐标
    vp_can = vp_can / vp_can(3);
    
    %检查所有的线段并计算inliers/outliers
     for k = 1 : num
        line_id_cur = line_ids(k);
        v1 = [lines(line_id_cur, 1), lines(line_id_cur, 3), 1];
        v2 = [lines(line_id_cur, 2), lines(line_id_cur, 4), 1];
        line = cross(v1, v2);
        line = line / line(3);
        error = abs(line*vp_can');
        if error <= eps
            cur_inliers = cur_inliers + 1;
            lines_inlier_cur = [lines_inlier_cur; lines(k, :)];
        else
            lines_outlier_cur = [lines_outlier_cur; lines(k, :)];
        end
     end
    
     %最适合
     if cur_inliers > max_inliers
        max_inliers = cur_inliers;
        vp = vp_can;
        lines_inlier_best = lines_inlier_cur;
        lines_outlier_best = lines_outlier_cur;
     end
end
        

end