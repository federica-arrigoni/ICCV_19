
function [missrate,known_rate,group_out,ind_known]=compute_missrate(group_out,labels_gt)
% Compute misclassification error and ratio of classified points
% Points labelled as unknown (if any) do not contribute to the error
% The range is [0,1]

% Classified points
ind_known=find(group_out~=0);

% Alignment 
group_out(ind_known)=bestMap(labels_gt(ind_known),group_out(ind_known));

% Compute misclassification error (only over classified points)
missrate = sum(labels_gt(ind_known) ~= group_out(ind_known)) / length(ind_known);
known_rate=length(ind_known)/length(labels_gt);


end
