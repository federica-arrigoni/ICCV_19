
% This script tests the "Mode" method on the dataset developed in:
% Robust motion segmentation from pairwise matches. Federica Arrigoni and
% Tomas Pajdla. ICCV 2019.  

clc,clear,close all
addpath(genpath('./'))

folder_path = './MY_DATASETS/';

dataset='FLOWERS';
%dataset='BAG';
%dataset='BEARS';
%dataset='PENCILS';
%dataset='PENGUIN';

img_path=[folder_path dataset '/'];

load([img_path 'data.mat'])
load([img_path 'gt_labels'])


%% Robust Motion Segmentation from Pairwise Matches - Mode

[group_mode]=segment_mode(labels_pairwise,pairwiseEst,dim,ncams,d);


%% Compute error

[missrate_mode,known_mode,group_mode,ind_mode]=compute_missrate(group_mode,labels_gt);

disp(['Missclassification error: ' num2str(missrate_mode*100) '%'])
disp(['Percentage of classified points: ' num2str(known_mode*100) '%'])


%% Visualize segmentation results

colors = lines(7);
colors=colors([4 3 2],:);

for i=1:ncams
    
    image_i=rgb2gray(imread(strcat(img_path,imnames(i).name))); % image i
    
    SIFT_i=SIFT{i};
    Xi=[SIFT_i.locs(:,1)';SIFT_i.locs(:,2)']'; % coordinates of points in image i
    ind_i=1+cumDim(i):cumDim(i)+dim(i); % indices of points in image i
    
    figure,
    imshow(image_i)
    hold on
    set(gca,'FontSize',22,'LineWidth',3)
    
    for p=1:length(ind_i) % points in image i
        if group_mode( ind_i(p) )~=0 % if the point has been classified...
            plot(Xi(p,1),Xi(p,2),'.','Color',colors(group_mode( ind_i(p) ),:),'MarkerSize',15)
        end
    end
    
end

%%


