% [Motion, Shape, T] = factorization(Wo, iterMax1, iterMax2)
% 
% This function computes the 3D object shape from missing and degenerate data. 
% 
% 
% Input arguments:
% 
%     Wo - The data matrix is defined as:
%                 Wo = [ u_1^1  ...  u_P^1
%                        v_1^1  ...  v_P^1
%                          .    .      .
%                          .     .     .
%                          .      .    .
%                        u_1^F  ...  u_P^F
%                        v_1^F  ...  v_P^F ]
%                 The missing data entries must be identified with NaN.
%                 
%     iterMax1 [Typical values: 3000-15000] - Maximum number of iterations for 
%                                             the missing data algorithm
%     iterMax2 [100-500] - Maximum number of iterations for Rigid Factorization
%     stopError1 [10^(-4) - 10^(-7)] - The missing data algorithm stops if the
%                                      error between two consecutive iterations
%                                      is lower than this threshold
%     stopError2 [10^(-1) - 10^(-3)] - Rigid Factorization algorithm stops if 
%                                      the error between two consecutive iterations 
%                                      is lower than this threshold
% 
% Note that Rigid Factorization is one step of the global missing data algorithm.
% 
% 
% Output arguments:
% 
%     Motion - The motion matrix stacks the camera matrices in all frames. The
%              camera matrix in frame f is a Stiefel matrix and is composed by lines 
%              2*f-1 and 2*f.
%     Shape - 3D object shape
%     T - Translation vector
%     
% 
% For details, see:
% 
%    Marques, M., and Costeira, J.. "Estimating 3D shape from degenerate sequences with missing data", 
%    Computer Vision and Image Understanding, 113(2):261-272, February 2009.
             
% NOTE: It has been optimized somewhat for the Reconstructing PASCAL VOC paper.

function [Motion, Shape, T, errors] = factorization_custom(Wo, iterMax1, iterMax2)


be_silent = true;


having_no_kp = find(all(isnan(Wo)'));
having_kps = setdiff(1:size(Wo,1), having_no_kp);
Wo(having_no_kp,:) = [];


WW = Wo;

M = not(isnan(Wo));
unknowns = sum(sum(not(M)));
Wo(find(M == 0)) = 0;
Wo = (sum(Wo,2)./sum(M,2)*ones(1,size(Wo,2))).*not(M) + Wo.*M;
W = Wo;

nImgs = floor(size(W,1)/2);
iter1 = 0;

diagError = [];

ind = find(sum(M,1) == 2*nImgs);

if length(ind) > 0

    T = Wo(:,ind(1));

else

    T = mean(W,2);
 
end

[o1,e,o2]=svd(W);
E=e(1:3,1:3);
O1=o1(:,1:3);
O2=o2(:,1:3);
R=O1*sqrtm(E);
S=sqrtm(E)*O2';

try,
while iter1 < iterMax1
    if ~be_silent
        t = tic();
        iter1
    end

    W_prev_J = W;

    W = W - T*ones(1,size(W,2));
    Woit = Wo - T*ones(1,size(W,2));

    Rprev = R;
    Wprev = W;

    iterAux = 0;

    Motion = zeros(size(Wo,1), 3);
    while iterAux < iterMax2
        %Motion = cell(nImgs,1);
        
        %parfor i = 1:nImgs                
        for i=1:nImgs
            range = 2*i-1:2*i;
            %%%% orig
            %A_f = projStiefel(R(2*i-1:2*i,:)');
            %%%% faster version
            [U,D,V] = svd(R(range,:)','econ');
            c = 0.5*(D(1,1) + D(2,2));
            A_f = c*U*V';

            %Motion{i} = A_f';
            Motion(range,:) = A_f';
        end
        %Motion = cell2mat(Motion);
        
        Shape = pinv(Motion)*W;
        R = W*pinv(Shape);

        Rprev = R;

        iterAux = iterAux + 1;        
    end

    W = Motion*Shape + T*ones(1,size(W,2));
    error1 = norm(W(M) - WW(M));

    change = norm(W(:) - W_prev_J(:));
    
    if ~be_silent
        fprintf('Projection error: %f, solution change: %f', error1, change);
    end

    W = Motion*Shape.*not(M) + Woit.*M + T*ones(1,size(Wo,2));        

    iter1 = iter1 + 1;

    if ~isempty(ind)
        T = Wo(:,ind(1));
    else
        T = mean(W,2);
        Motionret=Motion;
        Tret=T;
        Shaperet=Shape;
    end

    if ~be_silent
        toc(t)
    end
end
catch,
    Motion=Motioneret;
    T=Tret;
    Shape=Shaperet;
end

duh = Motion*Shape + T*ones(1,size(W,2));
counter = 1;
errors = zeros(((size(M,1)/2)),1);
for i=1:((size(M,1)/2))
    range = (2*(i-1) + 1):(2*i);
    this_duh = duh(range,:);
    this_duh = this_duh(M(range,:));
    this_WW = WW(range,:);
    this_WW = this_WW(M(range,:));
    % errors(counter) = norm(this_duh - this_WW);
    
    [sc1, scaling1] = scale_data(this_WW(1:2:end)', 'zscore');
    [sc2, scaling2] = scale_data(this_WW(2:2:end)', 'zscore');    
    this_sc1 = scale_data(this_duh(1:2:end)', 'zscore', scaling1);
    this_sc2 = scale_data(this_duh(2:2:end)', 'zscore', scaling2);

    % this is a custom error function, not the one that is optimized
    errors(counter) = norm([sc1; sc2] - [this_sc1; this_sc2]);
    counter = counter + 1;
end
newMotion = zeros(numel(having_kps) + numel(having_no_kp), 3);
newMotion(having_kps,:) = Motion;
newT = zeros(numel(having_kps ) + numel(having_no_kp),1);
newT(having_kps) = T;
Motion = newMotion;
T = newT;


function W = projStiefel(Wo)

[U,D,V] = svd(Wo,'econ');
c = mean(diag(D));

W = c*U*V';