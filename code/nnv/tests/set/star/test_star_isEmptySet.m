function test_star_isEmptySet()
    % TEST_STAR_ISEMPTYSET - Test Star.isEmptySet() method
    %
    % Tests that:
    %   1. A valid star set returns false (non-empty)
    %   2. An infeasible star set returns true (empty)

    % Create a non-empty star set
    center = [1; 1];
    V = [1 0 1; 0 1 1];
    P = ExamplePoly.randHrep('d', 3);
    S = Star([center V], P.A, P.b);

    % Check if non-empty
    b = S.isEmptySet;

    % Create an empty (infeasible) star set
    % Constraints: alpha <= 1 AND alpha >= 2 (impossible)
    S1 = Star([1 1], [1; -1], [1; -2]);
    b1 = S1.isEmptySet;

    % ASSERTIONS - verify correctness
    assert(b == 0, 'Star S should NOT be empty (isEmptySet should return 0)');
    assert(b1 == 1, 'Star S1 should be empty (isEmptySet should return 1)');

    % Save regression data
    data = struct();
    data.b_nonempty = b;
    data.b_empty = b1;
    data.S_V = S.V;
    data.S_C = S.C;
    data.S_d = S.d;
    data.S1_V = S1.V;
    data.S1_C = S1.C;
    data.S1_d = S1.d;
    save_test_data(data, 'test_star_isEmptySet', 'results', 'subdir', 'set');
end
