%% Visualize set creation using bias field, l_infinity and adjust contrast

% shared data
sliceSize = "96";
path2data = "../../data/UMCL/subjects/patient";
subject = "01";
c = 50; % around the middle of the volume image
% load 3d data
flair   = niftiread(path2data + subject+"/1/flair.nii");
flair = permute(flair, [3 2 1]);
flair = flair_normalization(flair);
mask    = niftiread(path2data + subject+"/1/mask.nii");
wm_mask = niftiread(path2data + subject+"/1/wm_mask.nii");

% Bias Field 
order = 3;
% coeff = [0.1, 0.25, 0.5];
% coeff_range = [0.00025, 0.0005, 0.001];
coeff = 0.5;
coeff_range = 0.001;

bbb = new_bias_field(img_norm, order, coeff);

[lb_bf, ub_bf] = generate_patches_bf(flair, mask, wm_mask, sliceSize, c, order, coeff, coeff_range); % generates all possible patches to analyze

mi = min(lb_bf, [], 'all');
ma = max(ub_bf,[], 'all');

f = figure('Position', get(0, 'Screensize'));
subplot(1,4,1)
imshow(squeeze(flair(c,:,:)), [min(flair(c,:,:), [], 'all'),max(flair(c,:,:),[], 'all')])
colorbar
title("Flair")

subplot(1,4,2)
imshow(lb_bf, [mi,ma])
colorbar
title("LB = BiasField(flair,coeff-coeffRange)")

subplot(1,4,3)
imshow(ub_bf, [mi,ma])
colorbar
title("UB = BiasField(flair,coeff+coeffRange)")

im_diff = ub_bf-lb_bf;
subplot(1,4,4)
imshow(im_diff, [min(im_diff, [], 'all'),max(im_diff,[], 'all')])
colorbar
title("UB - LB")

exportgraphics(f,'BiasFieldSample.pdf','ContentType','vector')


%% Adjust contrast variables

% gamma = [0.5; 1; 2]; % lower and upper bound for typical ranges used for gamma
% gamma_range = [0.0025; 0.00375; 0.005];
gamma = 0.5;
gamma_range = 0.005;

[lb_ac, ub_ac] = generate_patches_ac(flair, mask, wm_mask, sliceSize, c, gamma, gamma_range); % generates all possible patches to analyze

mi = min(lb_ac, [], 'all');
ma = max(ub_ac,[], 'all');

f = figure('Position', get(0, 'Screensize'));
subplot(1,4,1)
imshow(squeeze(flair(c,:,:)), [min(flair(c,:,:), [], 'all'),max(flair(c,:,:),[], 'all')])
colorbar
title("Flair")

subplot(1,4,2)
imshow(lb_ac, [mi,ma])
colorbar
title("LB = \gammaCorrection(flair,\gamma-\gammaRange)")

subplot(1,4,3)
imshow(ub_ac, [mi,ma])
colorbar
title("UB = \gammaCorrection(flair,\gamma+\gammaRange)")

im_diff = lb_ac-ub_ac;
subplot(1,4,4)
imshow(im_diff, [min(im_diff, [], 'all'),max(im_diff,[], 'all')])
colorbar
title("UB - LB")

exportgraphics(f,'gammaCorrectionSample.pdf','ContentType','vector')


%% L_inf variables
epsilon = 0.004;
nPix = 10; % percetnage of pixels to pertrub

[lb_l, ub_l] = generate_patches_linf(flair, mask, wm_mask, sliceSize, c, epsilon, nPix);

mi = min(lb_l, [], 'all');
ma = max(ub_l,[], 'all');

f = figure('Position', get(0, 'Screensize'));
subplot(1,4,1)
imshow(squeeze(flair(c,:,:)), [min(flair(c,:,:), [], 'all'),max(flair(c,:,:),[], 'all')])
colorbar
title("Flair")

subplot(1,4,2)
imshow(lb_l, [mi,ma])
colorbar
title("LB (Flair - \epsilon)")

subplot(1,4,3)
imshow(ub_l, [mi,ma])
colorbar
title("UB (Flair + \epsilon)")
im_diff = ub_l-lb_l;

subplot(1,4,4)
imshow(im_diff, [min(im_diff, [], 'all'),max(im_diff,[], 'all')])
colorbar
title("UB - LB")

exportgraphics(f,'randomNoiseSample.pdf','ContentType','vector')



%% Bias field functions
% Generate starting points for slices
function [flair_lb,flair_ub] = generate_patches_bf(flair, mask, wm_mask, sZ, c, order, coeffs, cRange)

    flair_slice   = squeeze(flair(c,:,:));

    % Apply transformation
    [flair_lb, flair_ub] = BiasField(flair_slice, order, coeffs, cRange);

end

% Get transformation bounds
function [img_lb, img_ub] = BiasField(img, order, coeffs, cRange)
% Get the code from torchIO to generate the BiasField
% https://torchio.readthedocs.io/_modules/torchio/transforms/augmentation/intensity/random_bias_field.html#RandomBiasField

    % Transform variables
    % coeffs = str2double(coeffs);
    % order  = str2double(order);
    % cRange = str2double(cRange);
    
    bField1 = generate_bias_field(img, coeffs-cRange, order);
    img1 = img .* bField1;

    bField2 = generate_bias_field(img, coeffs+cRange, order);
    img2 = img .* bField2;

    % get max and min value for every pixel given the biasField applied
    % interval range for every pixel is given by min and max values for
    % that pixel in images img1 and img2
    img_lb = min(img1,img2);
    img_ub = max(img1,img2);

    % Define input set as ImageStar
    % img_diff = img_ub - img_lb;
    % V(:,:,:,1) = img_lb; % assume lb is center of set (instead of img)
    % V(:,:,:,2) = img_diff ; % basis vectors
    % C = [1; -1]; % constraints
    % d = [1; -1];
    % I = ImageStar(V, C, d, 0, 1); % input set

end

% Generate bias field for transformation
function bias_field = generate_bias_field(img, coeffs, order)
% def generate_bias_field(
%     data: TypeData,
%     order: int,
%     coefficients: TypeData,
% ) -> np.ndarray:
%     # Create the bias field map using a linear combination of polynomial
%     # functions and the coefficients previously sampled
%     shape = np.array(data.shape[1:])  # first axis is channels
%     half_shape = shape / 2
% 
%     ranges = [np.arange(-n, n) + 0.5 for n in half_shape]
% 
%     bias_field = np.zeros(shape)
%     meshes = np.asarray(np.meshgrid(*ranges))
% 
%     for mesh in meshes:
%         mesh_max = mesh.max()
%         if mesh_max > 0:
%             mesh /= mesh_max
%     x_mesh, y_mesh, z_mesh = meshes
% 
%     i = 0
%     for x_order in range(order + 1):
%         for y_order in range(order + 1 - x_order):
%             for z_order in range(order + 1 - (x_order + y_order)):
%                 coefficient = coefficients[i]
%                 new_map = (
%                     coefficient
%                     * x_mesh**x_order
%                     * y_mesh**y_order
%                     * z_mesh**z_order
%                 )
%                 bias_field += np.transpose(new_map, (1, 0, 2))  # why?
%                 i += 1
%     bias_field = np.exp(bias_field).astype(np.float32)
%     return bias_field

    % Create the bias field map using a linear combination of polynomial
    % functions and the coefficients previously sampled
    shape = size(img); 
    % shape = shape(2:end); % first axis is channels (not necessary as it is a greyscale image -> 1 channel, already removed dimension)
    half_shape = shape ./ 2;
    ranges = {};
    
    for n = half_shape
        ranges{end+1} = (-n:1:(n-1)) + 0.5;
    end
    
    bias_field = zeros(shape);
    
    ndim = length(shape);
    meshes = zeros([ndim, shape(1), shape(2)]);
    for k=1:ndim
        meshes(k,:,:) = meshgrid(ranges{k});
    end
    
    for i = 1:size(meshes,1)
        mesh = meshes(i,:,:);
        mesh_max = max(mesh, [], 'all');
        if mesh_max > 0
            mesh = mesh./mesh_max;
            meshes(i,:,:) = mesh;
        end
    end
    
    x_mesh = meshes(1,:,:);
    y_mesh = meshes(2,:,:);
    
    % i = 0; % initialize counter (for coeffs)
    cf = coeffs;
    for x_order = 0:order
        for y_order = 0:(order-x_order)
            % cf = coeffs(i); % will this always be within limits? Why this?
            new_map = (cf .* x_mesh.^x_order .* y_mesh.^y_order);
            new_map = squeeze(permute(new_map, [2 1 3]));
            bias_field = bias_field + new_map; % why?
            % i = i + 1;
        end
    end
    
    % And that's it
    bias_field = exp(bias_field);

end

%% Adjust contrast functions

% Generate starting points for slices
function [flair_lb, flair_ub] = generate_patches_ac(flair, mask, wm_mask, sZ, c, gamma, gRange)

    flair_slice   = squeeze(flair(c,:,:));

    % Apply transformation
    [flair_lb, flair_ub] = AdjustContrast(flair_slice, gamma, gRange);

end

% Adjust contrast Perturbation
function [img1, img2] = AdjustContrast(img, gamma, gamma_range)
% Python Code from Project-MONAI
    % """
    % Apply the transform to `img`.
    % gamma: gamma value to adjust the contrast as function.
    % """
    % img = convert_to_tensor(img, track_meta=get_track_meta())
    % gamma = gamma if gamma is not None else self.gamma
    % 
    % if self.invert_image:
    %     img = -img
    % 
    % if self.retain_stats:
    %     mn = img.mean()
    %     sd = img.std()
    % 
    % epsilon = 1e-7
    % img_min = img.min()
    % img_range = img.max() - img_min
    % ret: NdarrayOrTensor = ((img - img_min) / float(img_range + epsilon)) ** gamma * img_range + img_min
    % 
    % if self.retain_stats:
    %     # zero mean and normalize
    %     ret = ret - ret.mean()
    %     ret = ret / (ret.std() + 1e-8)
    %     # restore old mean and standard deviation
    %     ret = sd * ret + mn
    % 
    % if self.invert_image:
    %     ret = -ret
    % 
    % return ret

    % Our Implementation
    img_min = min(img, [], 'all');
    img_max = max(img, [], 'all');
    img_range = img_max-img_min;
    epsilon = 1e-7; % not sure why we need this, but to avoid dividing by zero probably ???
    % gamma = str2double(gamma);
    % gamma_range = str2double(gamma_range);
    % if gamma == 2 % upper range
    %     gamma = gamma - gamma_range;
    % end
    % The range for gamma is 0.5 to 2 as default for most code I have seen
    
    % This is the transformed image. 
    % img_trans = ((img-img_min)./(img_range+epsilon)).^gamma * img_range + img_min;

    % Do we create the set from here by assuiming gama is a set of values?
    % This may make the most sense...
    % I believe this would assign no range values for the background
    % 
    % Would this look something like...
    img1 = ((img-img_min)./(img_range+epsilon)).^(gamma-gamma_range) * img_range + img_min;
    img2 = ((img-img_min)./(img_range+epsilon)).^(gamma+gamma_range) * img_range + img_min;
    % IS = ImageStar(lb,ub);
    
    % However, this will not be computationally efficient (defintely not as
    % efficient as the bright/dark perturbations)
    
    % Would this work?

    % Define input set as ImageStar
    % img_diff = img1 - img2;
    % V(:,:,:,1) = img2; % assume lb is center of set (instead of img)
    % V(:,:,:,2) = img_diff ; % basis vectors
    % C = [1; -1]; % constraints
    % d = [1; -1];
    % I = ImageStar(V, C, d, 0, 1); % input set

    % Is this okay? Check point
    % [lb1, ub1] = I.estimateRanges;
    % lb_diff = lb1 - img2;
    % ub_diff = ub1 - img1;

end


% Generate starting points for slices
function [flair_lb, flair_ub] = generate_patches_linf(flair, mask, wm_mask, sZ, c, epsilon, nPix)

    flair_slice   = squeeze(flair(c,:,:));
    mask_slice    = squeeze(mask(c,:,:));
    wm_mask_slice = squeeze(wm_mask(c,:,:));

    % Apply transformation
    [flair_lb, flair_ub, idxs] = L_inf(flair_slice, wm_mask_slice, epsilon, nPix);


end

% L_inf on randm pixels of each 2D patch
function [lb, ub, idxs] = L_inf(img, wm_mask, epsilon, nPix)

    rng(0); % to replicate results

    % Our Implementation
    % epsilon = str2double(epsilon);
    % nPix = str2double(nPix);

    % First get all the pixels in image that contains wm
    idxs = find(wm_mask == 1);
    N = length(idxs); % how many wm_pixels?
    cN = floor(N * nPix/100); % these are the total pixels to modify

    id = randperm(N,cN); % choose cN pixels out of N for each patch
    idxs = idxs(id);

    % Get data ranges
    img_min = min(img, [], 'all');
    img_max = max(img, [], 'all');
    img_range = img_max-img_min;

    epsilon = epsilon *img_range; % resize epsilon to match the equivalence of color values argued for

    % Apply perturnbation to those pixels
    lb = img;
    lb(idxs) = lb(idxs) - epsilon;
    ub = img;
    ub(idxs) = ub(idxs) + epsilon;

end

function bias_field_image = new_bias_field(img, degree, coeff)
    rank = length(size(img));
    coeff_mat = zeros(repmat(degree + 1, 1, rank));
    % Create coordinate vectors
    coords = cell(1, rank); % Initialize a cell array to hold coordinate vectors
    for i = 1:rank
        coords{i} = linspace(-1.0, 1.0, size(img,i));
    end
    % Get coefficients
    [row, col] = find(tril(ones(degree + 1)));
    linear_indices = sub2ind(size(coeff_mat), row, col);
    coeff_mat(linear_indices) = coeff;
    % Evaluate leggendre polynomials on x,y
    bf = evaluate_2d_legendre_series(coords{1}, coords{2}, coeff_mat); % Assuming 2D for leggrid2d
    % Apply field to image
    bias_field_image = img .* exp(bf)';
end

function result = evaluate_2d_legendre_series(x, y, c)
% Evaluates a 2-D Legendre series on the Cartesian product of x and y.
%
% Args:
%   x: 1-D array of x-coordinates.
%   y: 1-D array of y-coordinates.
%   c: 2-D array of coefficients, where c(i,j) is the coefficient of 
%      L_{i-1}(x) * L_{j-1}(y). Note the indexing shift due to Matlab's
%      1-based indexing.
%
% Returns:
%   result: 2-D array of the same size as meshgrid(x, y), containing the
%           values of the 2-D Legendre series.

    [X, Y] = meshgrid(x, y);  % Create the Cartesian product grid

    n_coeff_x = size(c, 1); % Number of coefficients in x-direction
    n_coeff_y = size(c, 2); % Number of coefficients in y-direction
    
    sum_val = 0;
    for k = 1:n_coeff_y  % Iterate over Legendre polynomials in y
        for l = 1:n_coeff_x  % Iterate over Legendre polynomials in x
            % Evaluate the Legendre polynomials and multiply by the coefficient
            sum_val = sum_val + c(l, k) * legendreP(l - 1, X) .* legendreP(k - 1, Y);
        end
    end
    result = sum_val;

end