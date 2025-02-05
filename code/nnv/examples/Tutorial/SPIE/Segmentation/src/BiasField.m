function bias_field_image = BiasField(img, degree, coeff)
    
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