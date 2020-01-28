
function [] = plot_sift(image_i,image_j,Zij,SIFT_i,SIFT_j)
% show corresponding points in two images

rows=size(image_i,1);
columns=50;
image_white=255*ones(rows,columns);

% ind1 and ind2 encode the matches between image_i and image_j
[ind1,ind2]=find(Zij);

% xi,xj are the matching points for the pair (i,j)
Xi=[SIFT_i.locs(ind1,1)';SIFT_i.locs(ind1,2)'];
Xj=[SIFT_j.locs(ind2,1)';SIFT_j.locs(ind2,2)'];

figure
if size(image_i,3)~=1
    image_i = rgb2gray(image_i);
    image_j = rgb2gray(image_j);
end

iptsetpref('ImshowInitialMagnification','fit');
imshow(cat(2, image_i, image_white, image_j)) ;

hold on ;
prl = lines(7);
for p=1:size(Xi,2) 
    
    color=prl(6,:); 
    plot([Xi(1,p) Xj(1,p)+ size(image_i,2)+columns],[Xi(2,p) Xj(2,p)],'-o','Color',color,'LineWidth',0.5,'MarkerSize',3);
end

axis image off ;


end