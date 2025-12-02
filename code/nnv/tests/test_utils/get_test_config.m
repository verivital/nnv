function config = get_test_config()
    % GET_TEST_CONFIG - Global test configuration for NNV test suite
    %
    % Usage:
    %   config = get_test_config();
    %   if config.save_figures
    %       save_test_figure(fig, 'test_name', 'fig_name', 1);
    %   end
    %
    % Configuration Options:
    %   save_figures         - Save figures to results/tests/figures/ (default: true)
    %   close_figures        - Close figures after saving (default: true)
    %   save_regression_data - Save .mat files for regression (default: true)
    %   compare_baselines    - Compare against saved baselines (default: false)
    %   tolerance            - Numerical comparison tolerance (default: 1e-6)
    %   results_dir          - Root directory for test results
    %   baselines_dir        - Directory for baseline .mat files
    %
    % Environment Variable Overrides:
    %   Set NNV_TEST_COMPARE_BASELINES=1 to enable baseline comparison (for CI/CD)
    %   Set NNV_TEST_SAVE_FIGURES=0 to disable figure saving
    %
    % To customize behavior, modify this file or override in individual tests:
    %   config = get_test_config();
    %   config.save_figures = false;  % Override for this test

    config = struct();

    % Figure saving options
    config.save_figures = true;       % Save figures on every run
    config.close_figures = true;      % Close figures after saving

    % Regression data options
    config.save_regression_data = true;  % Save .mat files for regression
    config.tolerance = 1e-6;             % Numerical comparison tolerance

    % Baseline comparison options
    % Default: OFF for local development, can be enabled for CI/CD
    config.compare_baselines = false;    % Compare against saved baselines
    config.fail_on_regression = true;    % Fail test if regression detected

    % Directory configuration
    config.results_dir = fullfile(nnvroot(), 'results', 'tests');
    config.figures_dir = fullfile(config.results_dir, 'figures');
    config.data_dir = fullfile(config.results_dir, 'data');
    config.baselines_dir = fullfile(config.results_dir, 'baselines');

    % Check for environment variable overrides (useful for CI/CD)
    if ~isempty(getenv('NNV_TEST_COMPARE_BASELINES'))
        config.compare_baselines = str2double(getenv('NNV_TEST_COMPARE_BASELINES')) == 1;
    end
    if ~isempty(getenv('NNV_TEST_SAVE_FIGURES'))
        config.save_figures = str2double(getenv('NNV_TEST_SAVE_FIGURES')) == 1;
    end
    if ~isempty(getenv('NNV_TEST_FAIL_ON_REGRESSION'))
        config.fail_on_regression = str2double(getenv('NNV_TEST_FAIL_ON_REGRESSION')) == 1;
    end
end
