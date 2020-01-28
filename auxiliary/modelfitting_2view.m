

function group=modelfitting_2view(y1,y2,ngroups)
% performs motion segmentation in two images by fitting multiple
% fundamental matrices to corresponding points with RPA
% Reference:
% Luca Magri and Andrea Fusiello. Robust multiple model fitting with
% preference analysis and low-rank approximation. BMVC 2015.

%% Parameters

sigma_gt=0.005; % default value

%% Set model

model = 'fundamental';
[ distFun, hpFun, fit_model, cardmss] = set_model( model );

%% normalization

npoints=size(y1,1);
y1=[y1';ones(1,npoints)];
y2=[y2';ones(1,npoints)];

[dat_img_1 T1] = normalise2dpts(y1);
[dat_img_2 T2] = normalise2dpts(y2);
X = [ dat_img_1 ; dat_img_2 ];

%% sampling hypotheses

% guided sampling
w = 0.5;
blk = 3*npoints;
S  = mssWeighted( X, 6*npoints, blk, 'cauchy', model, w, sigma_gt);
S=S(blk+1:end,:);

H = hpFun(X,S); %hypotheses

R = res( X, H, distFun ); disp('Residuals computed')

%% motion segmentation

group=model_fittingRPA(X,S,R,model,sigma_gt,ngroups);
    
end


function group=model_fittingRPA(X,S,R,model,sigma_gt,ngroups)


%% Preference Trick
P = prefMat(R, sigma_gt, 6); % preference matrix

K = exp(- (squareform(pdist(P,@tanimoto))).^2);  % similarity matrix

%% Robust PCA
try
    lambda = 1/sqrt(size(K,1));
    [K_rpca, E_hat, ~] = inexact_alm_rpca(K, lambda);
    
    
    %% symmetric matrix factorization
    
    %[Uinit, mekmeans]  = guess_init( K_rpca, ngroups , G);
    [Uinit]  = guess_init( K_rpca, ngroups);
    
    params.Hinit=Uinit; params.maxiter = 100000;
    [U, iter, obj] = symnmf_anls(K_rpca, ngroups,  params);
    indU = indMax(U);
    
    % segmentations obtained from snmf
    % NB: F is a segmentation
    F = seg_from_binaryU(U);
    softIndU= U;    softIndU(indU==0)=0; % mlsac like
    
    %% Model extraction
    
    niter = 1000;
    [Phi, Z] =  rinforzino(X, S, P, F, softIndU , model, sigma_gt, niter);
    
    [~,I] = max(Phi'*softIndU,[],1);
    mss = Z(I,:);
    
    
    %% refinement using robust statistic
    
    cost = 1.5271;
    group = segmentation( mss, X, model, U , sigma_gt, cost ,'nearest');
    
catch
    
    disp('Error in RPA: segmentation not performed')
    group=[];
    
end


end