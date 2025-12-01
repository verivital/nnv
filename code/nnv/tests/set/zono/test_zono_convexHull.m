function test_zono_convexHull()
    % TEST_ZONO_CONVEXHULL - Test Zono.convexHull() method
    %
    % Tests that:
    %   1. Convex hull produces valid output zonotope
    %   2. Hull contains all points from both input zonotopes

    % Create first zonotope
    c1 = [0; 0];
    V1 = [1 0 -1; 1 1 1];
    Z1 = Zono(c1, V1);

    % Create second zonotope
    c2 = [1; 1];
    V2 = [2 1 0; -1 1 0];
    Z2 = Zono(c2, V2);

    % Compute convex hull
    Z3 = Z2.convexHull(Z1);

    % ASSERTION 1: Result is valid
    assert(~isempty(Z3), 'Convex hull should be valid');
    assert(Z3.dim == 2, 'Convex hull should have dimension 2');

    % ASSERTION 2: Hull bounds should contain both input zonotopes
    tol = 1e-6;
    [lb1, ub1] = Z1.getBounds();
    [lb2, ub2] = Z2.getBounds();
    [lb3, ub3] = Z3.getBounds();

    expected_lb = min(lb1, lb2);
    expected_ub = max(ub1, ub2);

    assert(all(lb3 <= expected_lb + tol), ...
        'Convex hull lower bounds should be <= min of input lower bounds');
    assert(all(ub3 >= expected_ub - tol), ...
        'Convex hull upper bounds should be >= max of input upper bounds');

    % Create visualization
    fig = figure;
    Zono.plot(Z3);
    hold on;
    Zono.plot(Z2);
    hold on;
    Zono.plot(Z1);
    legend('Convex Hull Z3', 'Z2', 'Z1');
    title('Zonotope Convex Hull Test');

    save_test_figure(fig, 'test_zono_convexHull', 'convexHull', 1, 'subdir', 'set/zono');

    % Save regression data
    data = struct();
    data.Z1_c = Z1.c;
    data.Z1_V = Z1.V;
    data.Z2_c = Z2.c;
    data.Z2_V = Z2.V;
    data.Z3_c = Z3.c;
    data.hull_lb = lb3;
    data.hull_ub = ub3;
    save_test_data(data, 'test_zono_convexHull', 'results', 'subdir', 'set');
end
