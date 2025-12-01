function test_tutorial_figures()
    % TEST_TUTORIAL_FIGURES - Run tutorial examples and save all generated figures
    %
    % This test runs the tutorial examples that generate figures and saves
    % them to the test results directory for regression tracking.
    %
    % Figures saved:
    %   - MNIST input set examples (1 figure)
    %   - MNIST verification (2 figures)
    %   - MNIST FC verification (1 figure)
    %   - M2NIST segmentation (1 figure)
    %   - ACC verification (2 figures)
    %   - AEBS reachability (2 figures)
    %   - Inverted Pendulum reachability (1 figure)
    %   - Set representations (11 figures)
    %
    % Total: ~21 figures

    % Store original directory
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));

    % Base path for examples (nnvroot is nnv/, examples are in code/nnv/examples)
    examples_base = fullfile(nnvroot, 'code', 'nnv', 'examples');

    fig_count = 0;

    %% 1) MNIST: Input Set Examples
    try
        cd(fullfile(examples_base, 'Tutorial', 'NN', 'MNIST'));
        close all;
        input_set_examples;
        fig_count = fig_count + save_open_figures('tutorial_mnist_input_set', 'tutorial/mnist');
    catch ME
        warning('input_set_examples failed: %s', ME.message);
    end

    %% 2) MNIST: Verification
    try
        cd(fullfile(examples_base, 'Tutorial', 'NN', 'MNIST'));
        close all;
        verify;
        fig_count = fig_count + save_open_figures('tutorial_mnist_verify', 'tutorial/mnist');
    catch ME
        warning('%s', ME.message);
    end

    %% 3) MNIST: FC Verification
    try
        cd(fullfile(examples_base, 'Tutorial', 'NN', 'MNIST'));
        close all;
        verify_fc;
        fig_count = fig_count + save_open_figures('tutorial_mnist_verify_fc', 'tutorial/mnist');
    catch ME
        warning('%s', ME.message);
    end

    %% 4) Segmentation: M2NIST
    try
        cd(fullfile(examples_base, 'Tutorial', 'NN', 'Segmentation'));
        close all;
        verify_m2nist;
        fig_count = fig_count + save_open_figures('tutorial_segmentation_m2nist', 'tutorial/segmentation');
    catch ME
        warning('%s', ME.message);
    end

    %% 5) ACC: Verification
    try
        cd(fullfile(examples_base, 'Tutorial', 'NNCS', 'ACC', 'Verification'));
        close all;
        verify;
        fig_count = fig_count + save_open_figures('tutorial_acc_verify', 'tutorial/nncs');
    catch ME
        warning('%s', ME.message);
    end

    %% 6) AEBS: Reachability
    try
        cd(fullfile(examples_base, 'Tutorial', 'NNCS', 'AEBS'));
        close all;
        reach;
        fig_count = fig_count + save_open_figures('tutorial_aebs_reach', 'tutorial/nncs');
    catch ME
        warning('%s', ME.message);
    end

    %% 7) Inverted Pendulum: Reachability
    try
        cd(fullfile(examples_base, 'Tutorial', 'NNCS', 'InvertedPendulum'));
        close all;
        reach_invP;
        fig_count = fig_count + save_open_figures('tutorial_inverted_pendulum', 'tutorial/nncs');
    catch ME
        warning('%s', ME.message);
    end

    %% 8) Set Representations
    try
        cd(fullfile(examples_base, 'Tutorial', 'other'));
        close all;
        set_representations;
        fig_count = fig_count + save_open_figures('tutorial_set_representations', 'tutorial/sets');
    catch ME
        warning('%s', ME.message);
    end

    %% Summary
    fprintf('\n=== Tutorial Figure Test Summary ===\n');
    fprintf('Total figures saved: %d\n', fig_count);

    % Assertion to ensure we saved figures
    assert(fig_count >= 15, 'Expected at least 15 tutorial figures, got %d', fig_count);
end


function count = save_open_figures(base_name, subdir)
    % SAVE_OPEN_FIGURES - Save all currently open figures
    %
    % Inputs:
    %   base_name - Base name for saved figures
    %   subdir - Subdirectory under figures/tutorial/
    %
    % Returns:
    %   count - Number of figures saved

    count = 0;
    figs = findall(0, 'Type', 'figure');

    if isempty(figs)
        return;
    end

    for i = 1:length(figs)
        fig = figs(i);
        if ~isvalid(fig)
            continue;
        end

        % Generate figure name
        if length(figs) == 1
            fig_name = sprintf('%s', base_name);
        else
            fig_name = sprintf('%s_%d', base_name, i);
        end

        % Save using the test infrastructure
        save_test_figure(fig, base_name, fig_name, i, 'subdir', subdir);
        count = count + 1;
    end

    % Close all figures after saving
    close all;
end
