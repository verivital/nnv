function test_star_merge_stars()
    % TEST_STAR_MERGE_STARS - Test Star.merge_stars() method
    %
    % Tests that:
    %   1. Merging produces valid output stars
    %   2. Merged stars over-approximate the union of input stars

    % Create four star sets
    I1 = ExamplePoly.randVrep;
    V = [0 0; 1 0; 0 1];
    I1 = Star(V', I1.A, I1.b);

    I2 = ExamplePoly.randVrep;
    V = [1 1; 1 0; 0 1];
    I2 = Star(V', I2.A, I2.b);

    I3 = ExamplePoly.randVrep;
    V = [1 2; 1 0; 1 1];
    I3 = Star(V', I3.A, I3.b);

    I4 = ExamplePoly.randVrep;
    V = [1 1; 0 4; 3 1];
    I4 = Star(V', I4.A, I4.b);

    % Array of input stars
    I = [I1 I2 I3 I4];

    % Merge stars (merge into 2 groups)
    merge1 = Star.merge_stars(I, 2, 'single');

    % ASSERTION 1: Merged result is not empty
    assert(~isempty(merge1), 'merge_stars should return non-empty result');
    assert(length(merge1) == 2, 'Should have 2 merged star groups');

    % ASSERTION 2: Each merged star should be valid
    for i = 1:length(merge1)
        assert(~merge1(i).isEmptySet, sprintf('Merged star %d should not be empty', i));
    end

    % ASSERTION 3: Union of merged stars should over-approximate original stars
    tol = 1e-6;
    num_samples = 10;

    % Get bounds of all merged stars
    all_merged_lb = inf(I1.dim, 1);
    all_merged_ub = -inf(I1.dim, 1);
    for i = 1:length(merge1)
        [lb_i, ub_i] = merge1(i).getRanges();
        all_merged_lb = min(all_merged_lb, lb_i);
        all_merged_ub = max(all_merged_ub, ub_i);
    end

    % Verify original stars' samples are within merged bounds
    for j = 1:length(I)
        samples = I(j).sample(num_samples);
        for k = 1:size(samples, 2)
            s = samples(:, k);
            assert(all(s >= all_merged_lb - tol) && all(s <= all_merged_ub + tol), ...
                sprintf('Sample %d from star %d should be in merged bounds', k, j));
        end
    end

    % Create visualizations - merged stars
    fig1 = figure;
    Star.plot(merge1(1));
    title('Merged Star Group 1');
    save_test_figure(fig1, 'test_star_merge_stars', 'merged1', 1, 'subdir', 'set/star');

    fig2 = figure;
    Star.plot(merge1(2));
    title('Merged Star Group 2');
    save_test_figure(fig2, 'test_star_merge_stars', 'merged2', 2, 'subdir', 'set/star');

    % Create visualizations - original stars
    fig3 = figure;
    subplot(2, 2, 1); Star.plot(I1); title('Original I1');
    subplot(2, 2, 2); Star.plot(I2); title('Original I2');
    subplot(2, 2, 3); Star.plot(I3); title('Original I3');
    subplot(2, 2, 4); Star.plot(I4); title('Original I4');
    sgtitle('Original Stars Before Merge');
    save_test_figure(fig3, 'test_star_merge_stars', 'original', 3, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.num_input_stars = length(I);
    data.num_merged_stars = length(merge1);
    [data.merged1_lb, data.merged1_ub] = merge1(1).getRanges();
    [data.merged2_lb, data.merged2_ub] = merge1(2).getRanges();
    data.all_merged_lb = all_merged_lb;
    data.all_merged_ub = all_merged_ub;
    save_test_data(data, 'test_star_merge_stars', 'results', 'subdir', 'set');
end
