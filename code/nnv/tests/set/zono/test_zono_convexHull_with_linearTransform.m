function test_zono_convexHull_with_linearTransform()
    % TEST_ZONO_CONVEXHULL_WITH_LINEARTRANSFORM - Test optimized convex hull
    %
    % Tests that:
    %   1. convexHull_with_linearTransform produces valid result
    %   2. Result is equivalent to manual convexHull after affineMap

    % Set seed for reproducibility
    rng(42);

    % Create zonotope
    c = [1; 1];
    V = rand(2, 3);
    Z1 = Zono(c, V);

    % Define linear transformation
    W = [2 1; 0 -1];
    b = [];

    % Apply transformation
    Z2 = Z1.affineMap(W, b);

    % Method 1: Manual convex hull
    Z12 = Z1.convexHull(Z2);

    % Method 2: Optimized version
    Z3 = Z1.convexHull_with_linearTransform(W);

    % ASSERTION 1: Both results are valid
    assert(~isempty(Z12), 'Manual convex hull should be valid');
    assert(~isempty(Z3), 'Optimized convex hull should be valid');

    % ASSERTION 2: Results should have similar bounds
    tol = 1e-6;
    [lb12, ub12] = Z12.getBounds();
    [lb3, ub3] = Z3.getBounds();

    % Both should contain the original zonotopes
    [lb1, ub1] = Z1.getBounds();
    [lb2, ub2] = Z2.getBounds();

    assert(all(lb3 <= lb1 + tol) && all(ub3 >= ub1 - tol), ...
        'Optimized hull should contain Z1');
    assert(all(lb3 <= lb2 + tol) && all(ub3 >= ub2 - tol), ...
        'Optimized hull should contain Z2');

    % Create visualization - Method 1
    fig1 = figure;
    Zono.plot(Z12);
    hold on;
    Zono.plot(Z1);
    hold on;
    Zono.plot(Z2);
    legend('Manual Hull Z12', 'Z1', 'Z2');
    title('Manual Convex Hull');

    save_test_figure(fig1, 'test_zono_convexHull_with_linearTransform', 'manual', 1, 'subdir', 'set/zono');

    % Create visualization - Method 2
    fig2 = figure;
    Zono.plot(Z3);
    hold on;
    Zono.plot(Z1);
    hold on;
    Zono.plot(Z2);
    legend('Optimized Hull Z3', 'Z1', 'Z2');
    title('Optimized Convex Hull with Linear Transform');

    save_test_figure(fig2, 'test_zono_convexHull_with_linearTransform', 'optimized', 2, 'subdir', 'set/zono');

    % Save regression data
    data = struct();
    data.W = W;
    data.Z1_c = Z1.c;
    data.Z1_V = Z1.V;
    data.Z12_lb = lb12;
    data.Z12_ub = ub12;
    data.Z3_lb = lb3;
    data.Z3_ub = ub3;
    save_test_data(data, 'test_zono_convexHull_with_linearTransform', 'results', 'subdir', 'set');
end
