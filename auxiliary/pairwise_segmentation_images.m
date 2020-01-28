
function [labels_pairwise]=pairwise_segmentation_images(matches_pairwise,SIFT,d,n,A)
% Performs motion segmentation on image pairs with RPA (fitting multiple
% fundamental matrices to corresponding points)
%
% INPUT:
% matches_pairwise: pairwise matches (matches_pairwise{i,j}.ind1 and
% matches_pairwise{i,j}.ind2 contain the indices of corresponding points in
% image i and j) 
% SIFT: locations of image points
% n: number of images
% d: number of motions
% A: viewing graph
%
% OUTPUT:
% labels_pairwise: two-frame segmentations (labels_pairwise{i,j} contains
% the labels of corresponding points in images i and j)

if nargin<5
    A=ones(n); % viewing graph
end

labels_pairwise=cell(n);% initialization
for i=1:n
    
    SIFT_i=SIFT{i}; % points in image i
    
    for j=i+1:n
        
        if A(i,j)==1
            
            %%
            SIFT_j=SIFT{j}; % points in image j
            
            fprintf('\nPair of images %d %d: performing motion segmentation with RPA...\n',i,j) ;

            % ind1 and ind2 encode the matches between image_i and image_j
            ind1=matches_pairwise{i,j}.ind1;
            ind2=matches_pairwise{i,j}.ind2;
            
            % xi,xj are the matching points for the pair (i,j)
            Xi=[SIFT_i.locs(ind1,1)';SIFT_i.locs(ind1,2)'];
            Xj=[SIFT_j.locs(ind2,1)';SIFT_j.locs(ind2,2)'];
            
            if length(ind1)>8 % if there are (at least) 8 correspondences
                
                % perform motion segmentation in the pair (i,j) by fitting 
                % multiple fundamental matrices to correspoding points 
                groups=modelfitting_2view(Xi',Xj',d);
                
                % save labels
                labels_pairwise{i,j}=groups;
                labels_pairwise{j,i}=groups;
                
                fprintf('Done!\n') ;
            end
        end
    end
end


end




