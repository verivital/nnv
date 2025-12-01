function test_zono_orderReduction()
    % TEST_ZONO_ORDERREDUCTION - Test Zono.orderReduction_box() method
    %
    % Tests that:
    %   1. Order reduction produces valid zonotope
    %   2. Reduced zonotope over-approximates original

    % Create zonotope with high order
    c1 = [0; 0];
    V1 = [1 -1; 1 1; 0.5 1; -1.2 1];
    Z1 = Zono(c1, V1');

    % Get original order
    original_order = size(Z1.V, 2);

    % Apply order reduction (target order 3)
    target_order = 3;
    Z2 = Z1.orderReduction_box(target_order);

    % ASSERTION 1: Result is valid zonotope
    assert(~isempty(Z2), 'Reduced zonotope should be valid');
    assert(Z2.dim == Z1.dim, 'Reduced zonotope should have same dimension');

    % ASSERTION 2: Order should be reduced
    reduced_order = size(Z2.V, 2);
    assert(reduced_order <= target_order, ...
        sprintf('Reduced order (%d) should be <= target order (%d)', reduced_order, target_order));

    % ASSERTION 3: Reduced zonotope should over-approximate original
    tol = 1e-6;
    [lb1, ub1] = Z1.getBounds();
    [lb2, ub2] = Z2.getBounds();

    assert(all(lb2 <= lb1 + tol), ...
        'Reduced zonotope lower bounds should be <= original lower bounds');
    assert(all(ub2 >= ub1 - tol), ...
        'Reduced zonotope upper bounds should be >= original upper bounds');

    % Create visualization
    fig = figure;
    Zono.plot(Z2);
    hold on;
    Zono.plot(Z1);
    legend(sprintf('Reduced Z2 (order %d)', reduced_order), ...
           sprintf('Original Z1 (order %d)', original_order));
    title('Zonotope Order Reduction Test');

    save_test_figure(fig, 'test_zono_orderReduction', 'orderReduction', 1, 'subdir', 'set/zono');

    % Save regression data
    data = struct();
    data.Z1_c = Z1.c;
    data.Z1_V = Z1.V;
    data.Z2_c = Z2.c;
    data.Z2_V = Z2.V;
    data.original_order = original_order;
    data.reduced_order = reduced_order;
    data.Z1_lb = lb1;
    data.Z1_ub = ub1;
    data.Z2_lb = lb2;
    data.Z2_ub = ub2;
    save_test_data(data, 'test_zono_orderReduction', 'results', 'subdir', 'set');
end
