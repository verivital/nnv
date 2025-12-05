function test_zono_getSupInfinityNorm()
    % TEST_ZONO_GETSUPINFINITYNORM - Test Zono.getSupInfinityNorm() method
    %
    % Tests that:
    %   1. Supremum infinity norm is computed correctly
    %   2. Result is consistent with zonotope bounds

    % Create first zonotope
    c1 = [0; 0];
    V1 = [1 -1; 1 1];
    Z1 = Zono(c1, V1);
    r1 = Z1.getSupInfinityNorm();

    % Create second zonotope
    c2 = [1; 1];
    V2 = [2 1; -1 1];
    Z2 = Zono(c2, V2);
    r2 = Z2.getSupInfinityNorm();

    % ASSERTION 1: Norm should be positive
    assert(r1 > 0, 'Supremum infinity norm should be positive for Z1');
    assert(r2 > 0, 'Supremum infinity norm should be positive for Z2');

    % ASSERTION 2: Norm should be consistent with bounds
    % sup ||x||_inf = max(|lb|, |ub|) over all dimensions
    [lb1, ub1] = Z1.getBounds();
    [lb2, ub2] = Z2.getBounds();

    tol = 1e-6;
    expected_r1 = max(max(abs(lb1)), max(abs(ub1)));
    expected_r2 = max(max(abs(lb2)), max(abs(ub2)));

    assert(abs(r1 - expected_r1) < tol, ...
        'Z1 infinity norm should match expected value');
    assert(abs(r2 - expected_r2) < tol, ...
        'Z2 infinity norm should match expected value');

    % Create visualization
    fig = figure;
    Zono.plot(Z1);
    hold on;
    Zono.plot(Z2);
    legend(sprintf('Z1 (||x||_{inf} = %.2f)', r1), ...
           sprintf('Z2 (||x||_{inf} = %.2f)', r2));
    title('Zonotope Supremum Infinity Norm Test');

    save_test_figure(fig, 'test_zono_getSupInfinityNorm', 'infNorm', 1, 'subdir', 'set/zono');

    % Save regression data
    data = struct();
    data.Z1_c = Z1.c;
    data.Z1_V = Z1.V;
    data.r1 = r1;
    data.Z2_c = Z2.c;
    data.Z2_V = Z2.V;
    data.r2 = r2;
    save_test_data(data, 'test_zono_getSupInfinityNorm', 'results', 'subdir', 'set');
end
