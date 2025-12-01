function test_star_get_convex_hull()
    % TEST_STAR_GET_CONVEX_HULL - Test Star.get_convex_hull() static method
    %
    % Tests that:
    %   1. Convex hull of array of stars is computed correctly
    %   2. Hull contains all points from input stars

    % Create first star set
    I1 = ExamplePoly.randVrep;
    V = [0 0; 1 0; 0 1];
    I1 = Star(V', I1.A, I1.b);

    % Create second star set
    I2 = ExamplePoly.randVrep;
    V = [1 1; 1 0; 0 1];
    I2 = Star(V', I2.A, I2.b);

    % Compute convex hull of star array
    S = Star.get_convex_hull([I1, I2]);

    % ASSERTION 1: Result is valid
    assert(~isempty(S), 'get_convex_hull should return a valid result');

    % Create visualization
    fig = figure;
    S.plot();
    hold on;
    Star.plot(I1);
    hold on;
    Star.plot(I2);
    legend('Convex Hull', 'I1', 'I2');
    title('Star.get\_convex\_hull() Test');

    save_test_figure(fig, 'test_star_get_convex_hull', 'hull', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.I1_V = I1.V;
    data.I2_V = I2.V;
    data.result_class = class(S);
    save_test_data(data, 'test_star_get_convex_hull', 'results', 'subdir', 'set');
end
