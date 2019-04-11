load outputSet.mat;
load ACASXU_run2a_2_8_batch_2000.mat;

tic; 

R = R1;
n = length(R); 

% normalized matrix and vector
normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);


% normalize output set
parfor i=1:n
    fprintf('\nNormalize output set %d', i);
    R(i) = R(i).affineMap(normalized_mat, normalized_vec);   
end

% output: [x1 = COC; x2 = Weak Left; x3 = Weak Right; x4 = Strong Left; x5 = Strong Right]
% unsafe region: COC >= 15.8 & Strong Right <= 15.09

unsafe_mat = [-1 0 0 0 0; 0 0 0 0 1];
unsafe_vec = [-15.8; 15.09];

safe = 0;
US = []; % unsafe output set
parfor i=1:n
    S = R(i).intersectHalfSpace(unsafe_mat, unsafe_vec);
    if isempty(S)
        fprintf('\nThe %d^th output set does not reaches the unsafe region COC >= 15.8 & StrongLeft <= 15.09', i);
    else
        fprintf('\nThe %d^th output set reaches the unsafe region COC >= 15.8 & StrongLeft <= 15.09', i);
        safe = safe + 1;
        US = [US S];
    end
end

if safe >= 1
    fprintf('\n P41 is violated on N2_8');
else
    fprintf('\n P41 holds on N2_8');
end



% ========================================================================
% Input Constraints
lb = [1500; -0.06; 3.1; 1000; 700];
ub = [1800; 0.06; 3.14; 1200; 800];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end
B = Box(lb, ub);
I = B.toStar;
% =========================================================================


% construct a complete counter input set
n = length(US); % number of stars in counter input set
counterInputSet = [];
V = I.V; 
parfor i=1:n
   X = Star(V, US(i).C, US(i).d);
   counterInputSet = [counterInputSet X];
end


counterInputSet_constructionTime = toc;
save counterInputSet_constructionTime.mat counterInputSet_constructionTime;
save counterInputSet.mat counterInputSet;