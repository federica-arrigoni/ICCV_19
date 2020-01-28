function Q = permutation_procrustes( X, Y)
% Permutation Procrustes Analisys: it solves the problem Y=X*Q where Q is a
% permutation matrix using the Munkres algorithm. If X=I it finds the
% closest permutation. 

C = wthresh(-(X'*Y),'s', 0.0065);
Q = v2p(munkres(C),size(Y,2));
Q = spones(Q.*C); 
    
end


function [ P ] = v2p( v,c )
% transform vector into permutation matrix
n =  length(v);
if nargin < 2
    c = n;
end
i = (1:n);
i(v==0)=[];
v(v==0)=[];
P = sparse(i,v,1,n,c);

% check double stochastic
assert(all(sum(P,1)<=1) && all(sum(P,2)<=1))

end








