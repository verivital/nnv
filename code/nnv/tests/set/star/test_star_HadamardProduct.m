function test_star_HadamardProduct()
    % TEST_STAR_HADAMARDPRODUCT - Test Star.HadamardProduct() method
    %
    % Tests that:
    %   1. The product of two interval stars encloses the true product set,
    %      including the products attained at the interval endpoints
    %      (regression: the previous implementation returned [2, 2.5] for
    %      [1,2].*[1,2] and [39, 51] for [2,4].*[10,20], both disjoint
    %      from parts of the true product sets [1,4] and [20,80])
    %   2. The product is sound: for sampled x in S1, y in S2, the point
    %      x .* y is contained in S1.HadamardProduct(S2)
    %   3. Output star dimensions are consistent

    tol = 1e-6;

    %% 1a) Regression witness: product of [1,2] with itself
    % (two stars with identical constraints; previously the
    % shared-constraint branch returned ranges [2, 2.5])
    S = Star(1, 2);
    P = S.HadamardProduct(S);
    [lb, ub] = P.getRanges;
    assert(lb <= 1 + tol, 'Product of [1,2]x[1,2] must contain 1*1 = 1');
    assert(ub >= 4 - tol, 'Product of [1,2]x[1,2] must contain 2*2 = 4');

    %% 1b) Regression witness: x in [2,4], y in [10,20], scaled encodings
    % (different constraint matrices; previously the disjoint branch
    % returned ranges [39, 51], the true product set is [20, 80])
    S1 = Star([3 1],  [1; -1],     [1; 1]);      % x = 3 + a,   a in [-1,1]
    S2 = Star([15 5], [0.5; -0.5], [0.5; 0.5]);  % y = 15 + 5b, b in [-1,1]
    P2 = S1.HadamardProduct(S2);
    [lb2, ub2] = P2.getRanges;
    assert(lb2 <= 20 + tol, 'Product of [2,4]x[10,20] must contain 2*10 = 20');
    assert(ub2 >= 80 - tol, 'Product of [2,4]x[10,20] must contain 4*20 = 80');

    %% 2) Soundness on multi-dimensional stars with coupled predicates
    rng(0);
    V1 = [1 1 0; -2 1 1; 0.5 0 1];
    C1 = [1 0; -1 0; 0 1; 0 -1; 1 1];
    d1 = [1; 1; 1; 1; 1.5];
    A = Star(V1, C1, d1);

    V2 = [0 2 1; 3 0 1; -1 1 1];
    C2 = [1 0; -1 0; 0 1; 0 -1];
    d2 = [0.5; 0.5; 1; 1];
    B = Star(V2, C2, d2);

    PAB = A.HadamardProduct(B);
    assert(PAB.dim == A.dim, 'Product star must have the same dimension as the factors');
    assert(~PAB.isEmptySet, 'Product star should not be empty');

    [plb, pub] = PAB.getRanges;

    % bounds of the product star along random directions (catches
    % unsoundness that per-dimension ranges cannot)
    nDirs = 5;
    W = randn(nDirs, PAB.dim);
    wlb = zeros(nDirs, 1);
    wub = zeros(nDirs, 1);
    for k = 1:nDirs
        Pk = PAB.affineMap(W(k, :), []);
        [wlb(k), wub(k)] = Pk.getRanges;
    end

    % sample the factor predicates directly (rejection sampling, so no
    % polytope toolbox is needed) and check x .* y lands inside
    nChecked = 0;
    for trial = 1:200
        a = -1 + 2*rand(2, 1);
        if any(C1*a > d1)
            continue
        end
        b = -1 + 2*rand(2, 1);
        if any(C2*b > d2)
            continue
        end
        x = V1(:, 1) + V1(:, 2:end)*a;
        y = V2(:, 1) + V2(:, 2:end)*b;
        z = x .* y;
        nChecked = nChecked + 1;
        assert(all(z >= plb - tol) && all(z <= pub + tol), ...
            sprintf('Sampled product point (trial %d) must lie within the product star ranges', trial));
        wz = W*z;
        assert(all(wz >= wlb - tol) && all(wz <= wub + tol), ...
            sprintf('Sampled product point (trial %d) must lie within the product star projections', trial));
    end
    assert(nChecked > 50, 'Rejection sampling should accept a reasonable number of points');

    %% 3) Soundness across sign changes (factors straddling zero)
    S3 = Star(-1.5, 2);    % x in [-1.5, 2]
    S4 = Star(-3, 1);      % y in [-3, 1]
    P34 = S3.HadamardProduct(S4);
    [lb34, ub34] = P34.getRanges;
    % true product set is [-6, 4.5]
    assert(lb34 <= -6 + tol, 'Product of [-1.5,2]x[-3,1] must contain 2*(-3) = -6');
    assert(ub34 >= 4.5 - tol, 'Product of [-1.5,2]x[-3,1] must contain (-1.5)*(-3) = 4.5');
end
