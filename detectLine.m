%%���������ͼ���е�ֱ�ߣ��Ⱦ����Ե��Ϊ��һ������������׼��
%���룺
%grayIm: �������ĻҶ�ͼ��
%minLen: ��Ե����С���س���(����Ϊ0.025*�Խ���)
%�����
%lines��ֱ�ߵĲ������ṹΪ����nlines,[x1 x2 y1 y2 theta r]��

function lines = detectLine(grayIm,minLen)
%��ͼ������ÿһ�����ֵ�ݶ�
[dX, dY] = gradient(conv2(grayIm,fspecial('gaussian',7,15),'same'));

%ʹ��Canny�㷨���ұ�Ե
im_canny = edge(grayIm,'Canny');

%ȥ���߽��Ե
im_canny([1 2 end-1 end],:) = 0;
im_canny(:, [1 2 end-1 end]) = 0;

height = size(im_canny,1);
%���˵���Եͼ������ֵС�ڵ���0������ʵ���ǽ�ͼ��߽��������ȡ�����������ݴ˹����ݶ�,�����Ե����
ind = find(im_canny>0);
dX_ind = dX(ind);
dY_ind = dY(ind);
a_ind= atan(dY_ind./(dX_ind+1e-10));
num_dir = 8;
a_ind = ceil(mod(a_ind/pi*num_dir-0.5, num_dir));

% imshow(im_canny);

% ��ȡ������ߵ�����
for i =1:num_dir
    direction(i).ind = ind(find(a_ind == i));
end

%ɾ��̫С�ıߣ�������������ͬ�ı�һ��Ψһ��id
%��(�߶ȣ���ȣ�[�Ƕ�id])
lines = zeros(2000,6);

used = zeros(size(im_canny));

line_count = 0;
for k = 1:num_dir
    num_ind = 0;
    for m = (k-1):(k+1)
        num_ind = num_ind + sum(~used(direction(mod(m-1,num_dir)+1).ind));  %ͳ������3�������������
    end
    ind = zeros(num_ind,1);
    dir_im = zeros(size(im_canny));
    
    count =0;
    for m = (k-1):(k+1)
        m2 = mod(m-1,num_dir)+1;
        tind = direction(m2).ind(~used(direction(m2).ind));
        tmpcount = length(tind);
        ind(count+1:count+tmpcount) = tind;
        count = count + tmpcount;
    end
    
    dir_im(ind) = 1;
    [tmpL,num_edges] = bwlabel(dir_im,8);  %��8��ͨ������ͨ�򣬲����費ͬ��ͨ��ͬ��ǩ
    
    %��������ͨ��ֿ��洢��ע����ͨ����Կ���һ��ֱ�ߣ�
    edge_size = zeros(num_edges,1);
    edges = repmat(struct('ind',zeros(200,1)),num_edges,1);
    for i = 1:length(ind)
        id = tmpL(ind(i));
        edge_size(id) = edge_size(id)+1;
        edges(id).ind(edge_size(id)) = ind(i);
    end
    for i = 1:num_edges
        edges(i).ind = edges(i).ind(1:edge_size(i));
    end
    
    %��ȡ���ߵĶ˵�ͳ��ߵ�ͼ��
    for id = 1:num_edges
        if edge_size(id)>minLen
            y = mod(edges(id).ind-1,height)+1;      %ȷ��������ԭͼ�����������λ�ã�
            x = floor((edges(id).ind-1)/height)+1;  %ȷ��������ԭͼ�����������λ�ã�
            
            mean_x = mean(x);
            mean_y = mean(y);
            zmx = (x-mean_x);
            zmy = (y-mean_y);
            D = [sum(zmx.^2) sum(zmx.*zmy);sum(zmx.*zmy) sum(zmy.^2)];
            [v, lambda] = eig(D);
            theta = atan2(v(2,2),v(1,2));
            if lambda(1,1)>0
                conf = lambda(2,2)/lambda(1,1);
            else
                conf = 100000;
            end
            
            if conf>=400
                line_count = line_count+1;
                used(edges(id).ind) = 1;
                r = sqrt((max(x)-min(x))^2 + (max(y)-min(y))^2);
                x1 = mean_x - cos(theta)*r/2;
                x2 = mean_x + cos(theta)*r/2;
                y1 = mean_y - sin(theta)*r/2;
                y2 = mean_y + sin(theta)*r/2;
                
                r = mean_x*cos(theta)+mean_y*sin(theta);
                
                lines(line_count,1:6) = [x1 x2 y1 y2 theta r];
            end
        end
    end
end

lines = lines(1:line_count,:);
    
    
    
    
    
end


















