% test constructor for VolumeStar class

if ~is_github_actions() % out of memory error in github actions, need to update this with smaller images (create a small random one)

    load("img3d_data.mat"); % img_X and img_I
    
    % These are the shared variables we will use throughout this script
    
    disturbance = ones(size(img_I));
    
    % Constructor
    I = VolumeStar(img_I, -disturbance, disturbance);
    X = VolumeStar(img_X, -disturbance, disturbance);
    
    %% 1) Get range
    tic;
    [xmin, xmax] = I.getRange(1,1,1,1);
    toc;
    display(xmin);
    display(xmax);
    
    %% 2) AffineMap
    I2 = I.affineMap(1,1); % use scalars 
    
    % %% 3) Sampling
    % volumes = I.sample(10); % sample 10 images % TODO (takes too long, improve function)
    
    %% 4) Minkowski Sum
    I4 = I.MinkowskiSum(X); 
    
    %% 5) Concatenate
    I5 = I.concatenation(X);
    
    %% 6) Hadamard Product
    I6 = I.HadamardProduct(X);
    
    %% 7) isEmptySet
    res7 = I.isEmptySet;
    assert(res7 == 0);
    
    % %% 8) contains  %% TODO (out of memory, improve function)
    % res8 = I.contains(img_I);
    % assert(res8 == 1);
    % 
    % res8b = I.contains(img_X);
    % assert(res8 == 0);
    
    %% 9) estimateRanges
    [x1, x2] = I.estimateRanges;

end
