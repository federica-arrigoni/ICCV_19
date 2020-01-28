

function [pairwiseEst,Z] = matching_noransac(ncams,SIFT,dim)
% perform SIFT matching
% ncams: number of cameras
% dim(i) is the number of points in image i
% SIFT{i} contains SIFt locations/descriptors of image i
% pairwiseEst: pairwiseEst{i,j}.ind1 and pairwiseEst{i,j}.ind2 contain the
% indices of corresponding points in image i and image j
% Z is a block-matrix with permutation matrices representing matches

pairwiseEst = cell(ncams,ncams);
% pairwiseEst{i,j} contains the indices of the matches between image i and j

Z=sparse(sum(dim));
% Z contains relative permutations

cumDim = [0;cumsum(dim(1:end-1))];

for i = 1:ncams
    
    Zii=speye(dim(i));
    Z(1+cumDim(i):cumDim(i)+dim(i),1+cumDim(i):cumDim(i)+dim(i))=Zii;
    
    for j = i+1:ncams
        
        %% consider the pair (i,j)
        fprintf('\nPair of images %d %d\n',i,j) ;
        tic;
        [index,scores]=vl_ubcmatch(SIFT{i}.desc,SIFT{j}.desc,1.5); % match the pair (i,j)
        
        
        %% Match also the pair (j,i)
        [index_ji,scores_ji]=vl_ubcmatch(SIFT{j}.desc,SIFT{i}.desc,1.5); % match the pair (i,j)
        fprintf('Matched in %.3f s\n', toc) ;
        index_ji=fliplr(index_ji')';
        
        %% Remove non-symmetric matches
        ind_nomatch=[];
        for k=1:size(index,2)
            if ~ismember(index(:,k)',index_ji','rows')
                ind_nomatch=[ind_nomatch k];
            end
        end
        index(:,ind_nomatch)=[];
        
        %% Save matches and construct a permutation matrix
        
        % ind1 and ind2 encode the matches between image_i and image_j
        ind1=index(1,:);
        ind2=index(2,:);
        
%         % xi,xj are the matching points for the pair (i,j)
%         xi=[SIFT{i}.locs(ind1,1)';SIFT{i}.locs(ind1,2)'];
%         xj=[SIFT{j}.locs(ind2,1)';SIFT{j}.locs(ind2,2)'];
        
        fprintf('Number of matches: %d \n', length(ind1)) ;
        
        pairwiseEst{i,j}.ind1 = ind1;
        pairwiseEst{i,j}.ind2 = ind2;
        
        Zij=sparse(ind1,ind2,1,dim(i),dim(j));
        
        % chech that Zij is a permutation
        if ~isempty(find(sum(Zij)>1)) || ~isempty(find(sum(Zij')>1))
            fprintf('WARNING: no permutation for pair (%d,%d)\n',i,j) ;
        end
        
        Z(1+cumDim(i):cumDim(i)+dim(i),1+cumDim(j):cumDim(j)+dim(j))=Zij;
        Z(1+cumDim(j):cumDim(j)+dim(j),1+cumDim(i):cumDim(i)+dim(i))=Zij';
        
    end
end


end