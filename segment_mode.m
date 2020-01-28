

function [group]=segment_mode(labels_pairwise,matches_pairwise,dim,ncams,d)
% INPUT:
% labels_pairwise: two-frame segmentations (labels_pairwise{i,j} contains
% the labels of corresponding points in images i and j)
% matches_pairwise: pairwise matches (matches_pairwise{i,j}.ind1 and
% matches_pairwise{i,j}.ind2 contain the indices of corresponding points in
% image i and j) 
% dim: dim(i) is the number of points in image i
% ncams: number of images
% d: number of motions
%
% OUTPUT:
% group: labels of image points
%
% Author: Federica Arrigoni, 2019
% Reference: Robust motion segmentation from pairwise matches. Federica
% Arrigoni and Tomas Pajdla. ICCV 2019.

m=sum(dim); % total number of image points
group=zeros(m,1); % initialization
cumDim = [0;cumsum(dim(1:end-1))];


%% Count the number of classified points for each image pair

classified_points=zeros(ncams,ncams);
for i=1:ncams
    for j=i+1:ncams
        classified_points(i,j)=nnz(labels_pairwise{i,j} );
        classified_points(j,i)=classified_points(i,j);
        
        if isempty(labels_pairwise{j,i})
            labels_pairwise{j,i}=labels_pairwise{i,j};
        end
    end
end

% Construct a graph
A=(classified_points>=10);


%% Do segmentation for each image independently 

labels_absolute=cell(ncams,1);
for i=1:ncams
    
    estimates=[];
    
    % construct a graph where each node is a pair linked to image i
    nodes=find(A(i,:)); % number of pairs linked to image i
    n_nodes=length(nodes);
    
    % Compute pairwise permutations
    Q_pair=sparse(d*n_nodes,d*n_nodes);
    for k=1:n_nodes
        j=nodes(k);
        groups_ij=labels_pairwise{i,j}; % labels 
        Pk=labels2matrix(matches_pairwise,dim,i,j,groups_ij,d); % matrix representation
        
        for h=k+1:n_nodes
            j=nodes(h);
            groups_ij=labels_pairwise{i,j}; % labels
            Ph=labels2matrix(matches_pairwise,dim,i,j,groups_ij,d); % matrix representation
            
            % Solve an assignment problem
            Qhk=permutation_procrustes(Pk,Ph)';
            Q_pair(d*h-d+1:d*h,d*k-d+1:d*k)=Qhk;
            Q_pair(d*k-d+1:d*k,d*h-d+1:d*h)=Qhk';
        end
    end
    
    % Permutation synchronization
    Q_synch = pachauri_synch(Q_pair,d,n_nodes);
    
    % Update labels according to permutations
    for k=1:n_nodes
        j=nodes(k);
        groups_ij=labels_pairwise{i,j}; % labels of segmentation (i,j)
        [P_current]=labels2matrix(matches_pairwise,dim,i,j,groups_ij,d); % matrix representation
        
        Q=Q_synch(d*k-d+1:d*k,:); % optimal permutation
        P_current=P_current*Q;
        
        % Go back to labels
        labels_current=matrix2labels(matches_pairwise,i,j,groups_ij,P_current);
        [current_estimate]=labels_full(matches_pairwise,dim,i,j,labels_current); 
        
        % add one estimate for image i
        estimates=[estimates current_estimate];
    end

    
    %% Choose the most frequent label (mode) for each point
    
    n_observations=sum(~isnan(estimates),2); 
    
    % ignore all outliers
    estimates(estimates==0)=nan; 

    % compute the mode
    final_estimate=mode(estimates,2);
    final_estimate(isnan(final_estimate))=0;
    
    % count the number of entries that are equal to the mode: we require at
    % least two observations equal to the mode
    n_estimates=size(estimates,2);
    count=zeros(dim(i),1);
    for l=1:n_estimates
        count=count+( estimates(:,l)==final_estimate );
    end
    final_estimate( (count<=1) & (n_observations>=2) ) =0;
    
    % save labels of points in image i
    labels_absolute{i}=final_estimate;
    
end

%% Align all the segmentations

% Compute pairwise permutations
Q_pair=sparse(d*ncams,d*ncams);
for i=1:ncams
    
    Pi=a_labels2matrix(dim,i,labels_absolute,d); % matrix representation
    
    for j=i+1:ncams
        if A(i,j)==1
            
            Pij=labels2matrix(matches_pairwise,dim,i,j,labels_pairwise{i,j},d); % matrix representation
            Pji=labels2matrix(matches_pairwise,dim,j,i,labels_pairwise{i,j},d); % matrix representation
            Pj=a_labels2matrix(dim,j,labels_absolute,d); % matrix representation
            
            % Solve assignment problems
            Q = permutation_procrustes(Pij,Pi);           
            Q = Q'* permutation_procrustes(Pj,Pji)'; 
            
            Q_pair(d*i-d+1:d*i,d*j-d+1:d*j)=Q;
            Q_pair(d*j-d+1:d*j,d*i-d+1:d*i)=Q';
        end
    end
end

% Permutation synchronization
Q_synch = pachauri_synch(Q_pair,d,ncams);

% Update labels 
for j=1:ncams
    
    Q=Q_synch(d*j-d+1:d*j,:); % optimal permutation
    
    Pj=a_labels2matrix(dim,j,labels_absolute,d); % matrix representation
    Pj=Pj*Q; % Apply permutation
    
    % go back to labels
    final_labels=a_matrix2labels(Pj,dim,j,d);
    
    index_j=1+cumDim(j):cumDim(j)+dim(j); % points in image j
    group(index_j)=final_labels; % save labels
end


end

%%

function Q_synch = pachauri_synch(Q_pair,d,ncams)
% Perform permutation synchronization
% Reference: Pachauri et al. Solving the multi-way matching problem by
% permutation synchronization. NIPS 2013

% Compute the d leading eigenvectors
[U,Ds] = eigs(Q_pair,d,'lm');

% Rescale eigenvectors by corresponding eigenvalues
U=real(U)*sqrt(abs(Ds));

% Multiply by inverse of the first block
U1=U(1:d,:);
U=U*U1';

% Project onto permutations
Q_synch=zeros(d*ncams,d);
for i=1:ncams
    Q_synch(d*i-d+1:d*i,:)=permutation_procrustes(eye(d),U(d*i-d+1:d*i,:));
end

end


%%

function [P]=a_labels2matrix(dim,i,labels_absolute,d) 
% transform labels into a matrix representation (with reference to an
% absolute representation of segmentation)

labels=labels_absolute{i};
ind_inliers=find(labels~=0);
P=sparse(ind_inliers,labels(ind_inliers),1,dim(i),d); % image i

end

function labels=a_matrix2labels(P,dim,i,d)
% transform matrix representation into labels (with reference to an
% absolute representation of segmentation)

labels=zeros(dim(i),1); % image i
for k=1:d
    labels(find(P(:,k)))=k;
end

end

%%
function [P]=labels2matrix(matches_pairwise,dim,i,j,labels_reference,d)
% transform labels into a matrix representation (with reference to a local
% representation of segmentation)

% ind1 and ind2 encode the matches between image i and image j
if i<j
    ind1=matches_pairwise{i,j}.ind1;
    ind2=matches_pairwise{i,j}.ind2;
else
    ind1=matches_pairwise{j,i}.ind2;
    ind2=matches_pairwise{j,i}.ind1;
end

ind_inliers=find(labels_reference~=0);
P=sparse(ind1(ind_inliers),labels_reference(ind_inliers),1,dim(i),d); % image i

end

function labels=matrix2labels(matches_pairwise,i,j,groups_ij,P_current)
% transform labels into a matrix representation (with reference to a local
% representation of segmentation)

% ind1 and ind2 encode the matches between image i and image j
if i<j
    ind1=matches_pairwise{i,j}.ind1;
    ind2=matches_pairwise{i,j}.ind2;
else
    ind1=matches_pairwise{j,i}.ind2;
    ind2=matches_pairwise{j,i}.ind1;
end

labels=zeros(length(groups_ij),1);
for k=1:length(ind1)
    if ~isempty(find(P_current(ind1(k),:)))
        labels(k)=find(P_current(ind1(k),:));
    end
end

end

function [current_estimate]=labels_full(matches_pairwise,dim,i,j,labels_reference)
% extract labels of points in image i

% ind1 and ind2 encode the matches between image i and image j
if i<j
    ind1=matches_pairwise{i,j}.ind1;
    ind2=matches_pairwise{i,j}.ind2;   
else
    ind1=matches_pairwise{j,i}.ind2;
    ind2=matches_pairwise{j,i}.ind1;
end

current_estimate=nan( dim(i),1 );
current_estimate(ind1)=labels_reference;

end
