function result = evaluate_2d_legendre_series(x, y, c)

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