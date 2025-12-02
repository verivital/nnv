function [match, diff, details] = compare_regression_data(new_data, baseline_file, varargin)
    % COMPARE_REGRESSION_DATA - Compare new test output with baseline
    %
    % Usage:
    %   [match, diff] = compare_regression_data(new_data, 'baseline.mat')
    %   [match, diff] = compare_regression_data(new_data, 'baseline.mat', 'tolerance', 1e-8)
    %   [match, diff, details] = compare_regression_data(new_data, baseline_file)
    %
    % Inputs:
    %   new_data      - Struct containing new test results
    %   baseline_file - Path to baseline .mat file
    %
    % Optional Parameters:
    %   'tolerance' - Numerical tolerance for comparison (default: from config)
    %   'fields'    - Cell array of specific fields to compare (default: all common fields)
    %   'verbose'   - Print comparison details (default: false)
    %
    % Outputs:
    %   match   - Boolean indicating if all comparisons pass
    %   diff    - Struct with difference values for each compared field
    %   details - Cell array of {field, max_diff, passed} for each field
    %
    % Example:
    %   % Compare bounds from current run with baseline
    %   new_data = struct('lb', lb_computed, 'ub', ub_computed);
    %   [match, diff] = compare_regression_data(new_data, 'test_star_affineMap_bounds.mat');
    %   assert(match, 'Regression detected: bounds changed beyond tolerance');

    % Get default configuration
    config = get_test_config();

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'tolerance', config.tolerance, @isnumeric);
    addParameter(p, 'fields', {}, @iscell);
    addParameter(p, 'verbose', false, @islogical);
    parse(p, varargin{:});

    tolerance = p.Results.tolerance;

    % Initialize outputs
    match = true;
    diff = struct();
    details = {};

    % Check if baseline file exists
    if ~exist(baseline_file, 'file')
        warning('Baseline file not found: %s', baseline_file);
        match = false;
        diff.error = 'Baseline file not found';
        return;
    end

    % Load baseline data
    baseline = load(baseline_file);

    % Determine fields to compare
    if isempty(p.Results.fields)
        % Compare all common fields
        new_fields = fieldnames(new_data);
        baseline_fields = fieldnames(baseline);
        fields_to_compare = intersect(new_fields, baseline_fields);
    else
        fields_to_compare = p.Results.fields;
    end

    % Compare each field
    for i = 1:length(fields_to_compare)
        field = fields_to_compare{i};

        if ~isfield(new_data, field)
            warning('Field "%s" not found in new data', field);
            match = false;
            continue;
        end

        if ~isfield(baseline, field)
            warning('Field "%s" not found in baseline', field);
            match = false;
            continue;
        end

        new_val = new_data.(field);
        base_val = baseline.(field);

        % Check dimensions match
        if ~isequal(size(new_val), size(base_val))
            if p.Results.verbose
                fprintf('  %s: FAIL (dimension mismatch: %s vs %s)\n', ...
                    field, mat2str(size(new_val)), mat2str(size(base_val)));
            end
            diff.(field) = Inf;
            details{end+1, 1} = field;
            details{end, 2} = Inf;
            details{end, 3} = false;
            match = false;
            continue;
        end

        % Compute difference based on type
        if isnumeric(new_val) && isnumeric(base_val)
            field_diff = abs(new_val - base_val);
            max_diff = max(field_diff(:));
            passed = max_diff < tolerance;

            diff.(field) = max_diff;
            details{end+1, 1} = field;
            details{end, 2} = max_diff;
            details{end, 3} = passed;

            if ~passed
                match = false;
            end

            if p.Results.verbose
                if passed
                    fprintf('  %s: PASS (max diff: %.2e < %.2e)\n', field, max_diff, tolerance);
                else
                    fprintf('  %s: FAIL (max diff: %.2e >= %.2e)\n', field, max_diff, tolerance);
                end
            end

        elseif islogical(new_val) && islogical(base_val)
            passed = isequal(new_val, base_val);
            diff.(field) = ~passed;
            details{end+1, 1} = field;
            details{end, 2} = double(~passed);
            details{end, 3} = passed;

            if ~passed
                match = false;
            end

            if p.Results.verbose
                fprintf('  %s: %s (logical comparison)\n', field, ternary(passed, 'PASS', 'FAIL'));
            end

        else
            % For other types, use isequal
            passed = isequal(new_val, base_val);
            diff.(field) = ~passed;
            details{end+1, 1} = field;
            details{end, 2} = double(~passed);
            details{end, 3} = passed;

            if ~passed
                match = false;
            end

            if p.Results.verbose
                fprintf('  %s: %s (isequal comparison)\n', field, ternary(passed, 'PASS', 'FAIL'));
            end
        end
    end

    % Warn about fields only in new data
    if isempty(p.Results.fields)
        new_only = setdiff(fieldnames(new_data), baseline_fields);
        if ~isempty(new_only) && p.Results.verbose
            fprintf('  New fields (not compared): %s\n', strjoin(new_only, ', '));
        end
    end
end

function result = ternary(condition, if_true, if_false)
    if condition
        result = if_true;
    else
        result = if_false;
    end
end
