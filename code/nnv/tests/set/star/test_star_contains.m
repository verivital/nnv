function test_star_contains()
    % TEST_STAR_CONTAINS - Test Star.contains() method
    %
    % Tests that:
    %   1. Points inside the star return true (1)
    %   2. Points outside the star return false (0)

    % Create a star set from a random polytope
    I = ExamplePoly.randVrep;
    V = [0 0; 1 0; 0 1];
    S = Star(V', I.A, I.b);

    % Test point inside the star (should be contained)
    s1 = [-0.5; 1];

    % Test point outside the star (should not be contained)
    s2 = [1; 0.5];

    % Compute containment
    b1 = S.contains(s1);
    b2 = S.contains(s2);

    % ASSERTIONS - verify correctness
    assert(b1 == 1, 'Point s1 should be contained in star S');
    assert(b2 == 0, 'Point s2 should NOT be contained in star S');

    % Create visualization
    fig = figure;
    Star.plot(S);
    hold on;
    plot(s1(1,:), s1(2,:), 'go', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Inside (s1)');
    hold on;
    plot(s2(1,:), s2(2,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Outside (s2)');
    legend('Star S', 'Inside (s1)', 'Outside (s2)');
    title('Star Containment Test');

    % Save figure and data
    save_test_figure(fig, 'test_star_contains', 'containment', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.s1 = s1;
    data.s2 = s2;
    data.b1 = b1;
    data.b2 = b2;
    data.V = S.V;
    data.C = S.C;
    data.d = S.d;
    save_test_data(data, 'test_star_contains', 'results', 'subdir', 'set');
end
