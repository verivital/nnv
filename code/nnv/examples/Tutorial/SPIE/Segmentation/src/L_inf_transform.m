function [lb, ub, idxs] = L_inf_transform(img, wm_mask, epsilon, nPix)

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