load outputSet.mat;
load ACASXU_run2a_2_8_batch_2000.mat;

R = outputSet;
n = length(R); 

normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);

% normalize output set

parfor i=1:n
    fprintf('\nNormalize output set %d', i);
    R(i) = R(i).affineMap(normalized_mat, normalized_vec);
end

fprintf('\nGetting output range ...');
output_range = Star.get_hypercube_hull(R);
save output_range.mat output_range;
fprintf('\nGetting output range is done!');
