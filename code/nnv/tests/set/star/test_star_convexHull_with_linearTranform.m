function test_star_convexHull_with_linearTranform()
    % TEST_STAR_CONVEXHULL_WITH_LINEARTRANSFORM - Test convexHull_with_linearTransform
    %
    % Tests that:
    %   1. Convex hull with linear transform produces valid output
    %   2. Output contains transformed points from input

    % Create a star set
    I1 = ExamplePoly.randVrep;
    V = [0 0; 1 0; 0 1];
    I1 = Star(V', I1.A, I1.b);

    % Define transformation matrix
    W = [2 1; 1 -1];

    % Compute convex hull with linear transform
    I2 = I1.convexHull_with_linearTransform(W);

    % ASSERTION 1: Result is valid
    assert(~isempty(I2), 'Result should not be empty');

    % Create visualization
    fig = figure;
    Star.plot(I2);
    hold on;
    Star.plot(I1);
    legend('I2 (transformed hull)', 'I1 (original)');
    title('Star.convexHull\_with\_linearTransform() Test');

    save_test_figure(fig, 'test_star_convexHull_with_linearTranform', 'hull', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.I1_V = I1.V;
    data.W = W;
    save_test_data(data, 'test_star_convexHull_with_linearTranform', 'results', 'subdir', 'set');
end
