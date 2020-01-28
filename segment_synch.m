
function [group_synch] = segment_synch(Z,d,tau,theta)
% INPUT:
% Z: matrix representation of two-frame segmentations
% d: number of motions
% tau: threshold on the entries of the eigenvectors
% theta: threshold on the ratio between the largest and second-largest
% entry in each row of the eigenvectors
%
% OUTPUT:
% group_synch: labels of image points
%
% Author: Federica Arrigoni, 2019
% Reference: Motion segmentation via synchronization. Federica Arrigoni and
% Tomas Pajdla. ICCV Workshops 2019.


% Compute the d leading eigenvectors
[U,D] = eigs(Z,d,'la');

% Rescale eigenvectors by corresponding eigenvalues
U=real(U)*sqrt(abs(D));

% Ensure that the column sums are non-negative
U=U*diag(sign(sum(U))); 

% Set to zero rows smaller than tau
ind=find(max(U,[],2)<tau);
U(ind,:)=0;

% Set to zero rows such that the ratio between the largest entry and
% second-largest entry is smaller than theta
U_ratio = ratio_test(U,theta); 

% Compute the labels based on maximums over rows
[values,group_synch]=max((U_ratio),[],2); 
group_synch(values==0)=0;

end


function U = ratio_test(U,th)

for k=1:size(U,1)
    row=sort(abs(U(k,:)),'descend');
    
    if row(1)/row(2) <th
        U(k,:)=0;
    end
    
end

end


