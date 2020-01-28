
function [Z]=labels2Z(labels_pairwise,matches_pairwise,dim,d,n,A)
% construct a binary matrix Z that represent multiple two-frame
% segmentations

if nargin<6
    A=ones(n);
end

cumDim = [0;cumsum(dim(1:end-1))];
m=sum(dim);

Z=sparse(m,m); % initialization

for i=1:n
    for j=i+1:n
        if A(i,j)==1
            
            groups=labels_pairwise{i,j};
            
            if ~isempty(groups)
                
                % ind1 and ind2 encode the matches between image_i and image_j
                ind1=matches_pairwise{i,j}.ind1;
                ind2=matches_pairwise{i,j}.ind2;
                
                % assign a label to ALL points in image i (zero means no label)
                groups_i=zeros(dim(i),1);
                groups_i(ind1)=groups;
                
                % assign a label to ALL points in image j (zero means no label)
                groups_j=zeros(dim(j),1);
                groups_j(ind2)=groups;
                
                % construct a binary matrix that encodes the segmentation
                Zhk=segment2matrix(groups_i,groups_j,d);
                
                % update the block-matrix containing all segmentations
                Z(1+cumDim(i):cumDim(i)+dim(i),1+cumDim(j):cumDim(j)+dim(j)) = Zhk;
                Z(1+cumDim(j):cumDim(j)+dim(j),1+cumDim(i):cumDim(i)+dim(i)) = Zhk';
                
            end
        end
    end
end

end



function P = segment2matrix(x,y,d)
% d = number of motions 
% x(i)=k if point i in the first image belongs to motion k
% y(i)=k if point i in the second image belongs to motion k
% x(i)=0 if point i does not belong to any motion in common with y

n=length(x);
m=length(y);

P=zeros(n,m);

for i=1:d
    P(x==i,y==i)=1;
end


end


