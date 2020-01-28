

function pairwiseEst=ZtoMatches(Z,dim,ncams)
% pairwiseEst: pairwiseEst{i,j}.ind1 and pairwiseEst{i,j}.ind2 contain the
% indices of corresponding points in image i and image j
% Z is a block-matrix with permutation matrices representing matches
% ncams: number of cameras
% dim(i) is the number of points in image i

cumDim = [0;cumsum(dim(1:end-1))];
pairwiseEst = cell(ncams,ncams);

for i = 1:ncams
    for j=i+1:ncams
        
        % Zij is a permutation matrix representing the matches between
        % images i and j
        Zij=Z(1+cumDim(i):cumDim(i)+dim(i),1+cumDim(j):cumDim(j)+dim(j));
        [ind1,ind2]=find(Zij);
        
        pairwiseEst{i,j}.ind1 = ind1';
        pairwiseEst{i,j}.ind2 = ind2';
        
    end
end


end