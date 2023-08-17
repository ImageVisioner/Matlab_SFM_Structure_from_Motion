%%��������ȡ ����������,ʵ����������Զƽ��

function [vps] = detectVP(img_rgb)
img_gray = rgb2gray(img_rgb);
[heigh, width] = size(img_gray);

%Ѱ��ֱ��
lines = detectLine(img_gray, 0.025*sqrt(heigh*heigh+width*width));
[num_lines,~] = size(lines);
thetas = lines(:, 5);
idx = kmeans(thetas, 3);
color = ['r', 'g', 'b'];
figure(1);
hold off;
imshow(img_rgb);
hold on;

% ��ͬ�� id ��ֱ�߷ŵ�һ��������
clusters = {[],[],[]};
for i = 1:num_lines
    cluster_id = idx(i);
    clusters{cluster_id} = [clusters{cluster_id}, i];
    line(lines(i, [1,2])',lines(i,[3,4])','color',color(cluster_id));
end
vps = zeros(3, 3);
for i = 1:size(clusters,2)
    [vp,~] = ransacForVP(lines, clusters, i);
    vps(i, 1:2) = vp(1:2);
    vps(i, 3) = i;
end

end