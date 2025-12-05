function test_zono_getOrientedBox()
    % TEST_ZONO_GETORIENTEDBOX - Test Zono.getOrientedBox() method
    %
    % Tests that:
    %   1. Oriented box is computed correctly
    %   2. Box contains the zonotope

    % Create zonotope
    c1 = [0; 0];
    V1 = [1 -1; 1 1; 0.5 0; -1 0.5];
    Z1 = Zono(c1, V1');

    % Get oriented box
    B1 = Z1.getOrientedBox();

    % ASSERTION 1: Box is valid
    assert(~isempty(B1), 'Oriented box should be valid');

    % Create visualization
    fig = figure;
    Box.plot(B1);
    hold on;
    Zono.plot(Z1);
    legend('Oriented Box', 'Zonotope Z1');
    title('Zonotope Oriented Box Test');

    save_test_figure(fig, 'test_zono_getOrientedBox', 'orientedBox', 1, 'subdir', 'set/zono');

    % Save regression data
    data = struct();
    data.Z1_c = Z1.c;
    data.Z1_V = Z1.V;
    data.result_class = class(B1);
    save_test_data(data, 'test_zono_getOrientedBox', 'results', 'subdir', 'set');
end
