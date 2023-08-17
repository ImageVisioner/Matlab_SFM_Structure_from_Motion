function [IMGS] = im_read(count)
filePath = '..\images\';
IMGS = cell(count,1);
for i = 1:count
    id = i + 20;
    fileName = strcat(filePath,'B',int2str(id),'.jpg');
    IMGS{i} = imread(fileName);
end
end