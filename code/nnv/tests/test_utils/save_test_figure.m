function save_test_figure(fig_handle, test_name, fig_name, fig_num, varargin)
    % SAVE_TEST_FIGURE - Save test figure to results/tests/figures/
    %
    % Usage:
    %   save_test_figure(fig, 'test_star_affineMap', 'affineMap', 1)
    %   save_test_figure(fig, 'test_star_affineMap', 'affineMap', 1, 'save', false)
    %   save_test_figure(fig, 'test_star_affineMap', 'affineMap', 1, 'close', false)
    %   save_test_figure(fig, 'test_star_affineMap', 'affineMap', 1, 'subdir', 'set/star')
    %
    % Inputs:
    %   fig_handle  - Figure handle (gcf if empty)
    %   test_name   - Name of the test file (e.g., 'test_star_affineMap')
    %   fig_name    - Descriptive name for this figure (e.g., 'affineMap')
    %   fig_num     - Figure number (for tests with multiple figures)
    %
    % Optional Parameters:
    %   'save'   - Whether to save the figure (default: from config)
    %   'close'  - Whether to close the figure after saving (default: from config)
    %   'subdir' - Subdirectory under figures/ (e.g., 'set/star')
    %   'format' - Image format: 'png', 'fig', 'eps' (default: 'png')
    %
    % Output:
    %   Saves figure to: results/tests/figures/[subdir/]{test_name}_{fig_name}_{fig_num}.png

    % Get default configuration
    config = get_test_config();

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'save', config.save_figures, @islogical);
    addParameter(p, 'close', config.close_figures, @islogical);
    addParameter(p, 'subdir', '', @ischar);
    addParameter(p, 'format', 'png', @ischar);
    parse(p, varargin{:});

    % Use gcf if no handle provided
    if isempty(fig_handle)
        fig_handle = gcf;
    end

    % Save the figure if requested
    if p.Results.save
        % Determine output directory
        if isempty(p.Results.subdir)
            output_dir = config.figures_dir;
        else
            output_dir = fullfile(config.figures_dir, p.Results.subdir);
        end

        % Create directory if it doesn't exist
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end

        % Build filename
        filename = sprintf('%s_%s_%d.%s', test_name, fig_name, fig_num, p.Results.format);
        filepath = fullfile(output_dir, filename);

        % Save based on format
        switch lower(p.Results.format)
            case 'png'
                saveas(fig_handle, filepath);
            case 'fig'
                savefig(fig_handle, filepath);
            case 'eps'
                saveas(fig_handle, filepath, 'epsc');
            otherwise
                saveas(fig_handle, filepath);
        end

        % Print confirmation
        fprintf('  Saved figure: %s\n', filename);
    end

    % Close the figure if requested
    if p.Results.close
        close(fig_handle);
    end
end
