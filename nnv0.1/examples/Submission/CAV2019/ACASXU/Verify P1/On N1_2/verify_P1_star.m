load outputSet.mat;
load ACASXU_run2a_1_2_batch_2000.mat;

tic; 

R = outputSet;
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
% verify safety: COC <= 1500 or x1 <= 1500

unsafe_mat = [-1 0 0 0 0];
unsafe_vec = [-1500];

safe = 0;
parfor i=1:n
    S = R(i).intersectHalfSpace(unsafe_mat, unsafe_vec);
    if isempty(S)
        fprintf('\nThe %d^th output set does not reaches the unsafe region, P1 holds', i);
    else
        fprintf('\nThe %d^th output set reaches the unsafe region, P1 is violated', i);
        safe = safe + 1;
    end
end

if safe >= 1
    fprintf('\n P1 is violated on N1_2');
else
    fprintf('\n P1 holds on N1_2');
end

safety_checking_time = toc;
save safety_checking_time.mat safety_checking_time;
