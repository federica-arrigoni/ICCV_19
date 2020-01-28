
% This script tests the "Synch" method on the dataset developed in:
% Motion segmentation via synchronization. Federica Arrigoni and
% Tomas Pajdla. Workshop on Autonomous Navigation in Unconstrained 
% Environments (AUTONUE). ICCV Workshops 2019.  

clc,clear,close all
addpath(genpath('./'))

folder_path = './MY_DATASETS/';

dataset='CUPS';
%dataset='PEN';
%dataset='POUCH';
%dataset='BISCUITS';
%dataset='FOOD';
%dataset='NEEDLECRAFT';
%dataset='TEA';

img_path=[folder_path dataset '/'];

load([img_path 'data.mat'])
load([img_path 'gt_labels'])


%% Motion Segmentation via Synchronization - Synch

tau=0.01;
theta=1.5;
[group_synch] = segment_synch(Z,d,tau,theta);


%% Compute error

[missrate_synch,known_synch,group_synch,ind_synch]=compute_missrate(group_synch,labels_gt);

disp(['Missclassification error: ' num2str(missrate_synch*100) '%'])
disp(['Percentage of classified points: ' num2str(known_synch*100) '%'])


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
        if group_synch( ind_i(p) )~=0 % if the point has been classified...
            plot(Xi(p,1),Xi(p,2),'.','Color',colors(group_synch( ind_i(p) ),:),'MarkerSize',15)
        end
    end
    
end

%%



