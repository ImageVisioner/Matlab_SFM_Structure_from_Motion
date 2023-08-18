% 函数：获取SURF特征点和特征描述子
% 输入图片应为RGB图片，若不是，注释掉第15-17行代码
function [descriptors, points] = getPointFeatures(IMGS)

im_count = size(IMGS,1);
descriptors = cell(im_count,1);
points = cell(im_count,1);

if(im_count<2)
    msgbox('请输入足够的图像!','错误！');
    return;
end

for i =1:im_count
    if(ndims(IMGS{i})>2)
        image = rgb2gray(IMGS{i});
    end
    point = detectSURFFeatures(image);
    [descriptors{i},points{i}] = extractFeatures(image,point);
end
end