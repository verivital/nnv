function [sound, failed_points] = verify_soundness(input_set, output_set, transform_fn, varargin)
    % VERIFY_SOUNDNESS - Verify that a set transformation is sound
    %
    % A transformation is sound if: for all x in input_set, transform_fn(x) in output_set
    %
    % Usage:
    %   [sound, failed] = verify_soundness(S_in, S_out, @(x) W*x + b)
    %   [sound, failed] = verify_soundness(S_in, S_out, @(x) W*x + b, 'num_samples', 100)
    %
    % Inputs:
    %   input_set    - Input set (Star, Zono, Box, ImageStar, etc.)
    %   output_set   - Output set from the transformation
    %   transform_fn - Function handle that applies the transformation
    %
    % Optional Parameters:
    %   'num_samples' - Number of random samples to test (default: 50)
    %   'tolerance'   - Tolerance for containment check (default: 1e-6)
    %   'verbose'     - Print progress (default: false)
    %
    % Outputs:
    %   sound         - Boolean: true if all sampled points pass
    %   failed_points - Nx2 cell array of {input_point, output_point} that failed
    %
    % Example:
    %   % Verify affine transformation
    %   W = rand(3, 2);
    %   b = rand(3, 1);
    %   S_out = S_in.affineMap(W, b);
    %   [sound, failed] = verify_soundness(S_in, S_out, @(x) W*x + b);
    %   assert(sound, 'Affine transformation is not sound');

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'num_samples', 50, @isnumeric);
    addParameter(p, 'tolerance', 1e-6, @isnumeric);
    addParameter(p, 'verbose', false, @islogical);
    parse(p, varargin{:});

    num_samples = p.Results.num_samples;
    tolerance = p.Results.tolerance;
    verbose = p.Results.verbose;

    % Initialize outputs
    sound = true;
    failed_points = {};

    % Sample points from input set
    try
        if isa(input_set, 'Star')
            samples = input_set.sample(num_samples);
        elseif isa(input_set, 'Zono')
            % Zono doesn't have sample method, use getBox and sample from box
            [lb, ub] = input_set.getBounds();
            samples = zeros(length(lb), num_samples);
            for i = 1:num_samples
                samples(:, i) = lb + (ub - lb) .* rand(length(lb), 1);
            end
        elseif isa(input_set, 'Box')
            samples = zeros(length(input_set.lb), num_samples);
            for i = 1:num_samples
                samples(:, i) = input_set.lb + (input_set.ub - input_set.lb) .* rand(length(input_set.lb), 1);
            end
        elseif isa(input_set, 'ImageStar')
            % For ImageStar, sample in flattened space then reshape
            num_pixels = numel(input_set.V(:,:,:,1));
            samples_flat = input_set.sample(num_samples);
            samples = samples_flat;
        else
            warning('Unsupported set type: %s', class(input_set));
            sound = false;
            return;
        end
    catch ME
        warning('Failed to sample from input set: %s', ME.message);
        sound = false;
        return;
    end

    % Test each sample
    num_failed = 0;
    for i = 1:size(samples, 2)
        x = samples(:, i);

        % Apply transformation
        try
            y = transform_fn(x);
        catch ME
            if verbose
                fprintf('  Sample %d: transform failed - %s\n', i, ME.message);
            end
            num_failed = num_failed + 1;
            failed_points{end+1, 1} = x;
            failed_points{end, 2} = [];
            continue;
        end

        % Check containment in output set
        try
            if isa(output_set, 'Star')
                contained = output_set.contains(y);
            elseif isa(output_set, 'Zono')
                % For Zono, check if point is within bounds (over-approximation)
                [lb, ub] = output_set.getBounds();
                contained = all(y >= lb - tolerance) && all(y <= ub + tolerance);
            elseif isa(output_set, 'Box')
                contained = all(y >= output_set.lb - tolerance) && all(y <= output_set.ub + tolerance);
            else
                [lb, ub] = output_set.getRanges();
                contained = all(y >= lb - tolerance) && all(y <= ub + tolerance);
            end
        catch ME
            if verbose
                fprintf('  Sample %d: containment check failed - %s\n', i, ME.message);
            end
            % If containment check fails, try bounds-based check
            try
                [lb, ub] = output_set.getRanges();
                contained = all(y >= lb - tolerance) && all(y <= ub + tolerance);
            catch
                contained = false;
            end
        end

        if ~contained
            num_failed = num_failed + 1;
            failed_points{end+1, 1} = x;
            failed_points{end, 2} = y;
            sound = false;

            if verbose
                fprintf('  Sample %d: FAILED containment check\n', i);
            end
        end
    end

    if verbose
        fprintf('  Soundness: %d/%d samples passed\n', num_samples - num_failed, num_samples);
    end
end
