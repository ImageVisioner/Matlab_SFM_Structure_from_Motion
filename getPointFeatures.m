% ��������ȡSURF�����������������
% ����ͼƬӦΪRGBͼƬ�������ǣ�ע�͵���15-17�д���
function [descriptors, points] = getPointFeatures(IMGS)

im_count = size(IMGS,1);
descriptors = cell(im_count,1);
points = cell(im_count,1);

if(im_count<2)
    msgbox('�������㹻��ͼ��!','����');
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