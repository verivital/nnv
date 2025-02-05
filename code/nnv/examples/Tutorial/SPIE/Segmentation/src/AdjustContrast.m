function img = AdjustContrast(img, gamma)
% Code from Project-MONAI

    % Our Implementation
    img_min = min(img, [], 'all');
    img_max = max(img, [], 'all');
    img_range = img_max-img_min;
    epsilon = 1e-7; % not sure why we need this, but to avoid dividing by zero probably ???
    
    % transform image
    img = ((img-img_min)./(img_range+epsilon)).^(gamma) * img_range + img_min;

end