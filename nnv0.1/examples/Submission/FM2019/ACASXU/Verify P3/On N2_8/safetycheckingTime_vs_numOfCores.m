
load outputSet.mat;
load ACASXU_run2a_4_9_batch_2000.mat;

numCores = [1, 10, 20, 30, 40];
N = length(numCores);
safetycheckingTime = zeros(N,1); % safety checking time

for i=1:N
    
    % set up parallel computing with number of cores (workers)
    if numCores(i) >= 1
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            poolobj = parpool('local', numCores(i)); 
        else
            if poolobj.NumWorkers ~= numCores(i)
                delete(poolobj); % delete the old poolobj
                poolobj = parpool('local', numCores(i)); % start the new one with new number of cores
            end                    
        end
    end 
    
    tic; 

    R = R1;
    n = length(R); 

    % normalized matrix and vector
    normalized_mat = range_for_scaling(6) * eye(5);
    normalized_vec = means_for_scaling(6) * ones(5,1);


    % normalize output set
    parfor j=1:n
        fprintf('\nNormalize output set %d', j);
        R(j) = R(j).affineMap(normalized_mat, normalized_vec);   
    end

    % output: [x1 = COC; x2 = Weak Left; x3 = Weak Right; x4 = Strong Left; x5 = Strong Right]
    % safety property: COC is not the minimal score
    % unsafe region: COC is the minimal score: x1 <= x2; x1 <= x3; x1 <= x4, x1
    % <= x5

    unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
    unsafe_vec = [0; 0; 0; 0];

    safe = 0;
    parfor j=1:n
        S = R(j).intersectHalfSpace(unsafe_mat, unsafe_vec);
        if isempty(S)
            fprintf('\nThe %d^th output set does not reaches the unsafe region, P4 holds', j);
        else
            fprintf('\nThe %d^th output set reaches the unsafe region, P4 is violated', j);
            safe = safe + 1;
        end
    end

    if safe >= 1
        fprintf('\n P4 is violated on N2_8');
    else
        fprintf('\n P4 holds on N2_8');
    end

    safetycheckingTime(i) = toc;
    
end

save safetycheckingTime_star.mat safetycheckingTime;


