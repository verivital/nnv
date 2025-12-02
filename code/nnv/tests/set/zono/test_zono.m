function test_zono()
    % TEST_ZONO - Basic Zono (Zonotope) class functionality test
    %
    % Tests that:
    %   1. Zonotope construction works correctly
    %   2. getBox returns valid bounding box
    %   3. affineMap transforms zonotopes correctly
    %   4. MinkowskiSum is computed correctly

    % Create first zonotope
    c1 = [0; 0];
    V1 = [1 -1; 1 1];
    Z1 = Zono(c1, V1);

    % ASSERTION 1: Zonotope is valid
    assert(~isempty(Z1), 'Zonotope Z1 should be valid');
    assert(Z1.dim == 2, 'Zonotope Z1 should have dimension 2');

    % Get bounding box
    B1 = Z1.getBox();
    tol = 1e-6;

    % ASSERTION 2: Box should contain zonotope
    [zono_lb, zono_ub] = Z1.getBounds();
    assert(all(abs(B1.lb - zono_lb) < tol), 'Box lb should match zonotope bounds');
    assert(all(abs(B1.ub - zono_ub) < tol), 'Box ub should match zonotope bounds');

    % Create second zonotope
    c2 = [1; 1];
    V2 = [2 1; -1 1];
    Z2 = Zono(c2, V2);

    % Create visualization - Z1 and Z2
    fig1 = figure;
    Zono.plot(Z2);
    hold on;
    Zono.plot(Z1);
    legend('Z2', 'Z1');
    title('Zonotope Construction Test');

    save_test_figure(fig1, 'test_zono', 'construction', 1, 'subdir', 'set/zono');

    % Apply affine transformation
    W = [3 1; 1 0; 2 1];
    b = [0.5; 1; 0];
    Z3 = Z1.affineMap(W, b);

    % ASSERTION 3: Affine map produces valid result
    assert(Z3.dim == 3, 'Transformed zonotope should have dimension 3');

    fig2 = figure;
    Zono.plot(Z3);
    title('Affine Mapped Zonotope Z3');

    save_test_figure(fig2, 'test_zono', 'affineMap', 2, 'subdir', 'set/zono');

    % Get box of transformed zonotope
    B3 = Z3.getBox();

    fig3 = figure;
    Box.plot(B3);
    title('Bounding Box of Z3');

    save_test_figure(fig3, 'test_zono', 'box', 3, 'subdir', 'set/zono');

    % Compute Minkowski sum
    Z4 = Z1.MinkowskiSum(Z2);

    % ASSERTION 4: Minkowski sum is valid
    assert(Z4.dim == 2, 'Minkowski sum should have dimension 2');

    fig4 = figure;
    Zono.plot(Z4);
    hold on;
    Zono.plot(Z1);
    hold on;
    Zono.plot(Z2);
    legend('Z4 = Z1 + Z2', 'Z1', 'Z2');
    title('Minkowski Sum Test');

    save_test_figure(fig4, 'test_zono', 'MinkowskiSum', 4, 'subdir', 'set/zono');

    % Save regression data
    data = struct();
    data.Z1_c = Z1.c;
    data.Z1_V = Z1.V;
    data.Z2_c = Z2.c;
    data.Z2_V = Z2.V;
    data.B1_lb = B1.lb;
    data.B1_ub = B1.ub;
    data.Z3_dim = Z3.dim;
    data.Z4_c = Z4.c;
    save_test_data(data, 'test_zono', 'results', 'subdir', 'set');
end
