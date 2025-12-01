function test_star_getZono()
    % TEST_STAR_GETZONO - Test Star.getZono() method
    %
    % Tests that:
    %   1. getZono() returns a valid zonotope
    %   2. Zonotope over-approximates the star (contains all star points)

    % Set seed for reproducibility
    rng(3);

    % Create a star set
    center = [1; 1];
    V = [1 0 1; 0 1 1];
    P = ExamplePoly.randHrep('d', 3);
    S = Star([center V], P.A, P.b);

    % Get zonotope over-approximation
    Z = S.getZono();

    % Get oriented box
    B = S.getOrientedBox();

    % ASSERTION 1: Zonotope is valid
    assert(~isempty(Z), 'getZono should return a valid zonotope');
    assert(isa(Z, 'Zono'), 'Result should be a Zono object');

    % ASSERTION 2: Zonotope should over-approximate star
    % All star samples should be within zonotope bounds
    tol = 1e-6;
    num_samples = 20;
    samples = S.sample(num_samples);
    [zono_lb, zono_ub] = Z.getBounds();

    for i = 1:size(samples, 2)
        s = samples(:, i);
        assert(all(s >= zono_lb - tol) && all(s <= zono_ub + tol), ...
            sprintf('Star sample %d should be within zonotope bounds', i));
    end

    % ASSERTION 3: Star bounds should be within or equal to zonotope bounds
    [star_lb, star_ub] = S.getRanges();
    assert(all(zono_lb <= star_lb + tol), ...
        'Zonotope lower bounds should be <= star lower bounds');
    assert(all(zono_ub >= star_ub - tol), ...
        'Zonotope upper bounds should be >= star upper bounds');

    % Create visualization
    fig = figure;
    Zono.plot(Z);
    hold on;
    Box.plot(B);
    hold on;
    Star.plot(S);
    legend('Zonotope (getZono)', 'Oriented Box', 'Star');
    title('Star to Zonotope Conversion');

    save_test_figure(fig, 'test_star_getZono', 'zono', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.S_V = S.V;
    data.Z_c = Z.c;
    data.Z_V = Z.V;
    data.zono_lb = zono_lb;
    data.zono_ub = zono_ub;
    data.star_lb = star_lb;
    data.star_ub = star_ub;
    save_test_data(data, 'test_star_getZono', 'results', 'subdir', 'set');
end
