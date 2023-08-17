%% SFM测试文件
% 本测试文件，采用自标定的方法获取相机内参K，相关函数是 estimatingK;
% 采用DLT（直接线性求解）的方法处理PnP问题，并是用纯matlab编程实现
% 本文件主要为作者对3D重建中的“稀疏重建”的个人理解用笔记，有些地方采用了比较暴力的解决方法；
% 但此文件写出了解决SFM问题的大致流程，应对读者有一定的帮助
% 主要参考：
% 书籍：计算机视觉中多视图几何
% 源码：SFMedu; Homographies-for-Plane-Detection-and-3D-Reconstruction.

% 以下为主要代码，所用算法均放到algorithms文件夹下
clear
clc
% 加载所用函数所在文件的路径
addpath('F:\sfm-basic-process-master\algorithms');

%读取图片并提取特征点和特征描述子
images = imread("F:\sfm-basic-process-master\images\B22.jpg")
[descriptors,points] = getPointFeatures(images);

% 由第一幅图自标定相机内参K（注：假设所有图片的内参数一致）
num_img = size(images,1);
K = estimatingK(imread('F:\sfm-basic-process-master\images\B21.jpg'));

% 获取匹配点与各匹配点的索引
matchedPoints = cell(num_img-1,2);
matchPointsIndex = cell(num_img-1,1);
for i = 1:num_img-1
    indexPairs = matchFeatures(descriptors{i},descriptors{i+1});
    matchPointsIndex{i} = indexPairs; 
    matchedPoint1 = points{i}(indexPairs(:,1),:);
    matchedPoint2 = points{i+1}(indexPairs(:,2),:);
    matchedPoints{i,1} = matchedPoint1;
    matchedPoints{i,2} = matchedPoint2;
end

% 由前两幅图像匹配点估计基本矩阵F、本质矩阵E，并从本质矩阵中求相机2（对应第二幅图片）的外参，最后利用三角化，求匹配点的三维坐标
% 此处假设世界坐标系与相机1（对应第一幅图片）重合
F = estimateFundamentalMatrix(matchedPoints{1,1},matchedPoints{1,2});
E = computeE(F,K,K);
Rt1 = computePose(E,K,matchedPoints{1,1}.Location(1,:),matchedPoints{1,2}.Location(1,:));
X1 = trangulate(K,Rt1,matchedPoints{1,1},matchedPoints{1,2});
disp('Rt1 = ');
disp(Rt1);

% 求解PnP问题，即由3D点和其对应的2D点计算相机（对应图片）外参
% 为什么会有求解PnP的问题，以下是个人理解：
% 上面图像间匹配是邻近图像的匹配（例：image1和image2匹配、image2和image3匹配），
% 而世界坐标系与相机1重合，若是利用上面前两幅求外参的方法，只能的到相机3相对于相机2的外参，
% 并且此外参的位移向量是经过一个未知比例缩放过的，不是准确的。
% 本代码求PnP问题时，采用的是DLT的方式，为了避开BA优化（后面解释为什么避开），代码采用暴力匹配搜索并对应的方法
% 大致流程是：
% 1.找到第 i 副图像与第 i-1 副图像的匹配点和匹配点对应的3D坐标；
% 2.找到第 i 副图像与第 i+1 副图像的匹配点，并筛选出已知3D坐标的匹配点；
% 3.用DLT求解第 i+1 个相机（对应第 i+1 副图像）的外参；
% 4.筛选出第 i 副图像与第 i+1 副图像的未知3D坐标匹配点求取3D坐标并加入到已知3D坐标中去

% 相机3的外参（注意：此外参是世界坐标系下的，既是相对于相机1的）
[Rt2,X2,indForRt3] = getRt2andX(matchPointsIndex,points,X1,Rt1,K);

% 相机4的外参（注意：此外参是世界坐标系下的，既是相对于相机1的）
[Rt3,X3,indForRt4] = getRt3andX(matchPointsIndex,points,X2,Rt2,K,indForRt3);

% 相机5的外参（注意：此外参是世界坐标系下的，既是相对于相机1的）
[Rt4,X4] = getRt4andX(matchPointsIndex,points,X3,Rt3,K,indForRt4);

%% BA优化：重投影误差优化，是一种解决PnP问题的非线性优化方法。
% 在一般情况下，增量式的SFM在每次增加一幅图像（一个视角/一个相机），都需要BA优化，这会使代码运行时间较长
% 本代码的目的是了解SFM的基本流程，故选择了一种线性且较暴力的方法解决PnP的问题，可以提高代码运行速度
% BA优化在C++中使用ceres库实现较为简便，易于理解。而在matlab中实现较为复杂，不适合初学者阅读。
% 若是对BA优化的matlab代码有兴趣，可以参考SFMedu中的相关函数



            
        
        
        
        
        
        
        
        
        
        
        
        
        
        

