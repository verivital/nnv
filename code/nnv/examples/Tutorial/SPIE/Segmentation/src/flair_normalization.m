function [flair] = flair_normalization(flair)
% Normalize the images before verifying/inference

peak = get_peak(flair);
flair = flair ./ peak;

end

% get peak values based on KDE (for normalization)
function peak = get_peak(vol)

    % Get non-zero elements of vol and cast to double
    temp = double(vol(vol ~= 0));
    
    % Calculate 99th percentile
    q = prctile(temp, 99);
    
    % Keep values less than or equal to the 99th percentile
    temp = temp(temp <= q);
    
    % Reshape temp into a column vector
    temp = reshape(temp, [], 1);
    
    % Bandwidth for the KDE
    bw = q / 80;

    % Perform Kernel Density Estimation (KDE)
    [f, xi] = ksdensity(temp, 'Bandwidth', bw, 'npoints', 80, 'Kernel', 'normal');
    
    % Convert to percentage
    x_mat = 100.0 * f';
    y_mat = xi';
    
    % Find local maxima (peaks) in the density estimation
    [pks, locs] = findpeaks(x_mat);
    
    % Store the heights and locations of the peaks
    heights = pks;
    peaks = y_mat(locs);

    % get peak value
    [~, hfidx] = max(heights, [], 'all');
    peak = peaks(hfidx);

end