function filepath = save_test_data(data, test_name, data_name, varargin)
    % SAVE_TEST_DATA - Save test workspace data to results/tests/data/
    %
    % Usage:
    %   save_test_data(data, 'test_star_affineMap', 'affineMap_results')
    %   save_test_data(data, 'test_star_affineMap', 'results', 'subdir', 'set/star')
    %
    % Inputs:
    %   data      - Struct containing data to save
    %   test_name - Name of the test file (e.g., 'test_star_affineMap')
    %   data_name - Descriptive name for this data (e.g., 'results')
    %
    % Optional Parameters:
    %   'save'   - Whether to save the data (default: from config)
    %   'subdir' - Subdirectory under data/ (e.g., 'set/star')
    %
    % Output:
    %   filepath - Full path to saved file (empty if not saved)
    %   Saves data to: results/tests/data/[subdir/]{test_name}_{data_name}.mat
    %
    % Example:
    %   % Save reachability results for regression
    %   data = struct();
    %   data.input_lb = S.V(:,1) - norm(S.V(:,2:end), 'fro');
    %   data.input_ub = S.V(:,1) + norm(S.V(:,2:end), 'fro');
    %   data.output_lb = R.lb;
    %   data.output_ub = R.ub;
    %   save_test_data(data, 'test_star_affineMap', 'affineMap_bounds');

    filepath = '';

    % Get default configuration
    config = get_test_config();

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'save', config.save_regression_data, @islogical);
    addParameter(p, 'subdir', '', @ischar);
    parse(p, varargin{:});

    % Save the data if requested
    if p.Results.save
        % Determine output directory
        if isempty(p.Results.subdir)
            output_dir = config.data_dir;
        else
            output_dir = fullfile(config.data_dir, p.Results.subdir);
        end

        % Create directory if it doesn't exist
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end

        % Build filename
        filename = sprintf('%s_%s.mat', test_name, data_name);
        filepath = fullfile(output_dir, filename);

        % Save the data struct
        save(filepath, '-struct', 'data');

        % Print confirmation
        fprintf('  Saved data: %s\n', filename);
    end
end
