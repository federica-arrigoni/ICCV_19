
% This is the demo for running the "Mode" method on raw images.
% Reference: 
% Robust motion segmentation from pairwise matches. Federica Arrigoni and
% Tomas Pajdla. ICCV 2019.  

clc,clear,close all 
addpath(genpath('./'))

folder_path = './MY_DATASETS/';

dataset='FLOWERS'; d=2; % number of motions
%dataset='BAG'; d=2; % number of motions
%dataset='BEARS'; d=3; % number of motions
%dataset='PENCILS'; d=2; % number of motions
%dataset='PENGUIN'; d=2; % number of motions

format_img='JPG';
img_path=[folder_path dataset '/'];

imnames = dir(strcat(img_path,'*.',format_img));
ncams=length(imnames);


%% Parameters

scale=1; % (no scaling)
%scale=0.5; % Rescale images to speed up SIFT

% delete features that have <= min_match matches
min_match=1;

%% Compute SIFT locations and descriptors for each image

SIFT = cell(1,ncams);
dim=zeros(ncams,1);

for i=1:ncams
    
    fprintf('\nComputing frames and descriptors: image %d \n',i);
    
    tic;
    im = imread(strcat(img_path,imnames(i).name)); % load the current image
    im = imresize(im,scale); % rescale the image to speed-up SIFT
    
    if size(im,3)==1
        im=single(im);
    else
        im=single(rgb2gray(im));
    end
    [frames1,descr1] = vl_sift(im) ; % computes SIFT locations and descriptors
    
    SIFT{i}.desc = descr1;
    SIFT{i}.locs = frames1(1:2,:)';
    SIFT{i}.locs = SIFT{i}.locs/scale;
    SIFT{i}.scale = frames1(3,:)';
    
    fprintf('%d descriptors extracted\n',size(SIFT{i}.locs,1));
    toc
    
    dim(i)=size(SIFT{i}.locs,1); % number of points in image i
    
end

cumDim = [0;cumsum(dim(1:end-1))];

%% Match all the image pairs

[pairwiseEst,Z_pairwise] = matching_noransac(ncams,SIFT,dim);

% delete features that have <= min_match matches
for i=1:ncams
    n_match=sum( Z_pairwise(1+cumDim(i):cumDim(i)+dim(i),:) ,2);
    ind_match=find(n_match<=min_match);
    
    Z_pairwise(cumDim(i)+ind_match,:)=[];
    Z_pairwise(:,cumDim(i)+ind_match)=[];
    dim(i)=dim(i)-length(ind_match);
    
    cumDim = [0;cumsum(dim(1:end-1))];
    
    SIFT{i}.desc(:,ind_match)=[];
    SIFT{i}.locs(ind_match,:)=[];
    SIFT{i}.scale(ind_match)=[];
end

m=size(Z_pairwise,1); % total number of image points
pairwiseEst=ZtoMatches(Z_pairwise,dim,ncams);

%% Plot matches for a single pair

i=1; j=3;

image_i=imread(strcat(img_path,imnames(i).name)); % left image
image_j=imread(strcat(img_path,imnames(j).name)); % right image
SIFT_i=SIFT{i};
SIFT_j=SIFT{j};

Zij=Z_pairwise(1+cumDim(i):cumDim(i)+dim(i),1+cumDim(j):cumDim(j)+dim(j));
plot_sift(image_i,image_j,Zij,SIFT_i,SIFT_j);

set(gca,'FontSize',22,'LineWidth',3)
title(['Pair (' num2str(i) ',' num2str(j) ')'],'FontWeight','Normal')


%% Perform pairwise segmentation

% A is the viewing graph 
% It says in which image pairs there are matches and there is motion
if strcmp(dataset,'PENGUIN')
    A=ones(ncams);
    i=1; j=2; A(i,j)=0; A(j,i)=0;
    i=3; j=4; A(i,j)=0; A(j,i)=0;
    i=5; j=6; A(i,j)=0; A(j,i)=0;
else
    A=ones(ncams);
end

[labels_pairwise]=pairwise_segmentation_images(pairwiseEst,SIFT,d,ncams,A);


%% Plot segmentation for an image pair

l=1; r=3;

image_l=imread(strcat(img_path,imnames(l).name));
image_r=imread(strcat(img_path,imnames(r).name));
if size(image_l,3)~=1
    image_l = rgb2gray(image_l);
    image_r = rgb2gray(image_r);
end
[rows,columns]=size(image_r);

SIFT_l=SIFT{l}; SIFT_r=SIFT{r};
Xl=[SIFT_l.locs(:,1)';SIFT_l.locs(:,2)']'; % points in image l
Xr=[SIFT_r.locs(:,1)';SIFT_r.locs(:,2)']'; % points in image r

Zij=Z_pairwise(1+cumDim(l):cumDim(l)+dim(l),1+cumDim(r):cumDim(r)+dim(r));
[ind1,ind2]=find(Zij);
Xl=Xl(ind1,:); Xr=Xr(ind2,:); % mathcing points in images l,r

ind_l=1+cumDim(l):cumDim(l)+dim(l); % indices of points in image i
ind_r=1+cumDim(r):cumDim(r)+dim(r); % indices of points in image i

colors = lines(7);

group_lr=labels_pairwise{l,r};

figure,
imshow(cat(2, image_l, image_r)) ;
hold on
set(gca,'FontSize',22,'LineWidth',3)
title(['Input - pair (' num2str(l) ',' num2str(r) ')'],'FontWeight','Normal')

for p=1:length(ind1)
    plot(Xl(p,1),Xl(p,2),'.','Color',colors(group_lr( p )+1,:),'MarkerSize',15)
    plot(Xr(p,1)+columns,Xr(p,2),'.','Color',colors(group_lr( p )+1,:),'MarkerSize',15)
end
% blu is for outliers


%% Motion Segmentation from Pairwise Matches - ICCV19

[group_out]=segment_mode(labels_pairwise,pairwiseEst,dim,ncams,d);

% Classified points
ind_known=find(group_out~=0);
known_rate=length(ind_known)/length(group_out);
disp(['Percentage of classified points: ' num2str(known_rate*100) '%'])

%% Motion segmentation via synchronization - AUTONUE19

% [Z]=labels2Z(labels_pairwise,pairwiseEst,dim,d,ncams);
% tau=0.01;
% theta=1.5;
% [group_out] = segment_synch(Z,d,tau,theta);
% 
% % Classified points
% ind_known=find(group_out~=0);
% known_rate=length(ind_known)/length(group_out);
% disp(['Percentage of classified points: ' num2str(known_rate*100) '%'])

%% Plot segmentation for each image

colors = lines(7);
colors=colors([4 3 2],:);

for i=1:ncams
    
    image_i=imread(strcat(img_path,imnames(i).name)); % image i
    if size(image_i,3)~=1
        image_i = rgb2gray(image_i);
    end
    
    SIFT_i=SIFT{i};
    Xi=[SIFT_i.locs(:,1)';SIFT_i.locs(:,2)']'; % coordinates of points in image i
    ind_i=1+cumDim(i):cumDim(i)+dim(i); % indices of points in image i
    
    figure,
    imshow(image_i)
    hold on
    set(gca,'FontSize',22,'LineWidth',3)
    
    for p=1:length(ind_i) % points in image i
        if group_out( ind_i(p) )~=0 % if the point has been classified
            plot(Xi(p,1),Xi(p,2),'.','Color',colors(group_out( ind_i(p) ),:),'MarkerSize',15)
        end
    end
    
    
end



