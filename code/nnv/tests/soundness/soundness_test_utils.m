classdef soundness_test_utils
    % SOUNDNESS_TEST_UTILS Helper functions for soundness verification tests
    % These utilities verify that computed reachable sets contain all actual
    % outputs for sampled inputs from the input set.

    methods (Static)

        function contained = verify_star_containment(S, point, tol)
            % VERIFY_STAR_CONTAINMENT Check if a point is contained in a Star set
            %
            % Input:
            %   S - Star set
            %   point - column vector to check
            %   tol - tolerance for numerical comparisons (default: 1e-6)
            %
            % Output:
            %   contained - true if point is in S, false otherwise

            if nargin < 3
                tol = 1e-6;
            end

            % Star representation: x = c + V*alpha where C*alpha <= d
            c = S.V(:, 1);
            V_basis = S.V(:, 2:end);

            n_pred = size(V_basis, 2);
            n_out = length(c);

            if n_pred == 0
                % No predicate variables - just check if point equals center
                contained = norm(point - c) < tol;
                return;
            end

            % Solve: find alpha such that V_basis * alpha = point - c
            %        subject to C * alpha <= d
            residual = point - c;

            % For under-determined systems (n_pred > n_out), use LP to find
            % a feasible alpha that satisfies both equality and inequality constraints
            if n_pred > n_out && ~isempty(S.C)
                % Solve feasibility LP
                f = zeros(n_pred, 1);  % No objective
                Aeq = V_basis;
                beq = residual;
                A = S.C;
                b = S.d;
                lb = S.predicate_lb;
                ub = S.predicate_ub;

                options = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');
                [alpha_solve, ~, exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, options);

                if exitflag == 1
                    reconstruction_error = norm(V_basis * alpha_solve - residual);
                    contained = reconstruction_error < tol;
                else
                    contained = false;
                end
            else
                % Standard case: solve linear system
                if rank(V_basis) < n_pred
                    alpha_solve = pinv(V_basis) * residual;
                else
                    alpha_solve = V_basis \ residual;
                end

                % Check reconstruction error
                reconstruction_error = norm(V_basis * alpha_solve - residual);
                if reconstruction_error > tol
                    contained = false;
                    return;
                end

                % Check constraint satisfaction
                if isempty(S.C)
                    contained = true;
                else
                    contained = all(S.C * alpha_solve <= S.d + tol);
                end
            end
        end

        function contained = verify_imagestar_containment(IS, image, tol)
            % VERIFY_IMAGESTAR_CONTAINMENT Check if an image is contained in an ImageStar
            %
            % Input:
            %   IS - ImageStar set
            %   image - image array to check (same dimensions as IS center)
            %   tol - tolerance (default: 1e-6)
            %
            % Output:
            %   contained - true if image is in IS, false otherwise

            if nargin < 3
                tol = 1e-6;
            end

            % Direct check without toStar conversion
            n_pred = IS.numPred;

            if n_pred == 0
                % No predicate variables - check if image equals center
                contained = norm(image(:) - IS.V(:,:,:,1)) < tol;
                return;
            end

            % Flatten everything for linear algebra
            center = reshape(IS.V(:,:,:,1), [], 1);
            point = reshape(image, [], 1);

            V_basis = zeros(length(center), n_pred);
            for k = 1:n_pred
                V_basis(:, k) = reshape(IS.V(:,:,:,k+1), [], 1);
            end

            % Solve: V_basis * alpha = point - center
            residual = point - center;
            n_out = length(center);

            % For under-determined systems (n_pred > n_out), use LP to find
            % a feasible alpha that satisfies both equality and inequality constraints
            if n_pred > n_out && ~isempty(IS.C)
                % Solve feasibility LP:
                % min 0 subject to:
                %   V_basis * alpha = residual (equality)
                %   C * alpha <= d (inequality)
                %   pred_lb <= alpha <= pred_ub (bounds)

                f = zeros(n_pred, 1);  % No objective

                % Equality constraints: V_basis * alpha = residual
                Aeq = V_basis;
                beq = residual;

                % Inequality constraints: C * alpha <= d
                A = IS.C;
                b = IS.d;

                % Bounds
                lb = IS.pred_lb;
                ub = IS.pred_ub;

                % Solve LP
                options = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');
                [alpha_solve, ~, exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, options);

                if exitflag == 1
                    % Found feasible solution - verify reconstruction
                    reconstruction_error = norm(V_basis * alpha_solve - residual);
                    contained = reconstruction_error < tol;
                else
                    % No feasible solution exists
                    contained = false;
                end
            else
                % Standard case: solve linear system
                if rank(V_basis) < n_pred
                    alpha_solve = pinv(V_basis) * residual;
                else
                    alpha_solve = V_basis \ residual;
                end

                % Check reconstruction error
                reconstruction_error = norm(V_basis * alpha_solve - residual);
                if reconstruction_error > tol
                    contained = false;
                    return;
                end

                % Check constraint satisfaction
                if isempty(IS.C)
                    contained = true;
                else
                    contained = all(IS.C * alpha_solve <= IS.d + tol);
                end
            end
        end

        function samples = sample_star(S, n)
            % SAMPLE_STAR Generate n random samples from a Star set
            %
            % Input:
            %   S - Star set
            %   n - number of samples to generate
            %
            % Output:
            %   samples - matrix where each column is a sample point

            c = S.V(:, 1);
            V_basis = S.V(:, 2:end);
            n_pred = size(V_basis, 2);

            if n_pred == 0
                samples = repmat(c, 1, n);
                return;
            end

            samples = zeros(length(c), n);

            for i = 1:n
                % Generate random alpha within bounds
                alpha = S.predicate_lb + (S.predicate_ub - S.predicate_lb) .* rand(n_pred, 1);

                % Check constraint satisfaction and retry if needed
                max_tries = 100;
                tries = 0;
                while ~isempty(S.C) && any(S.C * alpha > S.d) && tries < max_tries
                    alpha = S.predicate_lb + (S.predicate_ub - S.predicate_lb) .* rand(n_pred, 1);
                    tries = tries + 1;
                end

                samples(:, i) = c + V_basis * alpha;
            end
        end

        function samples = sample_imagestar(IS, n)
            % SAMPLE_IMAGESTAR Generate n random samples from an ImageStar set
            %
            % Input:
            %   IS - ImageStar set
            %   n - number of samples to generate
            %
            % Output:
            %   samples - cell array of sample images

            img_size = size(IS.V);
            img_size = img_size(1:end-1);  % Remove predicate dimension

            n_pred = IS.numPred;
            samples = cell(1, n);

            for i = 1:n
                if n_pred == 0
                    samples{i} = IS.V(:,:,:,1);
                else
                    % Generate random alpha within bounds
                    alpha = IS.pred_lb + (IS.pred_ub - IS.pred_lb) .* rand(n_pred, 1);

                    % Check constraint satisfaction and retry if needed
                    max_tries = 100;
                    tries = 0;
                    while ~isempty(IS.C) && any(IS.C * alpha > IS.d) && tries < max_tries
                        alpha = IS.pred_lb + (IS.pred_ub - IS.pred_lb) .* rand(n_pred, 1);
                        tries = tries + 1;
                    end

                    % Compute image: center + sum(alpha_i * basis_i)
                    img = IS.V(:,:,:,1);
                    for j = 1:n_pred
                        img = img + alpha(j) * IS.V(:,:,:,j+1);
                    end
                    samples{i} = img;
                end
            end
        end

        function corners = get_predicate_corners(n_pred, pred_lb, pred_ub)
            % GET_PREDICATE_CORNERS Get all corner points of the predicate hypercube
            %
            % Input:
            %   n_pred - number of predicate variables
            %   pred_lb - lower bounds (n_pred x 1)
            %   pred_ub - upper bounds (n_pred x 1)
            %
            % Output:
            %   corners - matrix where each column is a corner (n_pred x 2^n_pred)

            if n_pred == 0
                corners = [];
                return;
            end

            % Limit to avoid explosion
            if n_pred > 10
                warning('Too many predicate variables for corner sampling, using random');
                corners = pred_lb + (pred_ub - pred_lb) .* rand(n_pred, min(1024, 2^n_pred));
                return;
            end

            n_corners = 2^n_pred;
            corners = zeros(n_pred, n_corners);

            for i = 1:n_corners
                bits = bitget(i-1, 1:n_pred);
                for j = 1:n_pred
                    if bits(j) == 0
                        corners(j, i) = pred_lb(j);
                    else
                        corners(j, i) = pred_ub(j);
                    end
                end
            end
        end

        function [passed, msg] = verify_layer_soundness_star(layer, input_star, method, n_samples, tol)
            % VERIFY_LAYER_SOUNDNESS_STAR Verify soundness of layer reach for Star input
            %
            % Input:
            %   layer - NNV layer object with evaluate and reach methods
            %   input_star - Star input set
            %   method - reach method (e.g., 'approx-star', 'exact-star')
            %   n_samples - number of random samples to test
            %   tol - tolerance
            %
            % Output:
            %   passed - true if all samples contained in reach result
            %   msg - error message if failed

            if nargin < 4
                n_samples = 50;
            end
            if nargin < 5
                tol = 1e-5;
            end

            % Compute reachable set
            output_sets = layer.reach(input_star, method);
            % Handle different return types: single set, array, or cell array
            if iscell(output_sets)
                % Already a cell array
            elseif length(output_sets) > 1
                % Array of Stars - convert to cell array
                temp = cell(1, length(output_sets));
                for idx = 1:length(output_sets)
                    temp{idx} = output_sets(idx);
                end
                output_sets = temp;
            else
                % Single Star
                output_sets = {output_sets};
            end

            passed = true;
            msg = '';

            % Test corner cases
            corners = soundness_test_utils.get_predicate_corners(...
                input_star.nVar, input_star.predicate_lb, input_star.predicate_ub);

            for i = 1:size(corners, 2)
                alpha = corners(:, i);

                % Check constraint satisfaction
                if ~isempty(input_star.C) && any(input_star.C * alpha > input_star.d + tol)
                    continue;  % Skip infeasible corners
                end

                % Compute concrete input and output
                input_concrete = input_star.V(:,1) + input_star.V(:,2:end) * alpha;
                output_concrete = layer.evaluate(input_concrete);

                % Check containment in any output set
                contained = false;
                for j = 1:length(output_sets)
                    if soundness_test_utils.verify_star_containment(output_sets{j}, output_concrete, tol)
                        contained = true;
                        break;
                    end
                end

                if ~contained
                    passed = false;
                    msg = sprintf('Corner %d not contained in output', i);
                    return;
                end
            end

            % Test random samples
            samples = soundness_test_utils.sample_star(input_star, n_samples);

            for i = 1:n_samples
                input_concrete = samples(:, i);
                output_concrete = layer.evaluate(input_concrete);

                contained = false;
                for j = 1:length(output_sets)
                    if soundness_test_utils.verify_star_containment(output_sets{j}, output_concrete, tol)
                        contained = true;
                        break;
                    end
                end

                if ~contained
                    passed = false;
                    msg = sprintf('Random sample %d not contained in output', i);
                    return;
                end
            end
        end

        function [passed, msg] = verify_layer_soundness_imagestar(layer, input_is, method, n_samples, tol)
            % VERIFY_LAYER_SOUNDNESS_IMAGESTAR Verify soundness of layer reach for ImageStar
            %
            % Input:
            %   layer - NNV layer object with evaluate and reach methods
            %   input_is - ImageStar input set
            %   method - reach method (e.g., 'approx-star', 'exact-star')
            %   n_samples - number of random samples to test
            %   tol - tolerance
            %
            % Output:
            %   passed - true if all samples contained in reach result
            %   msg - error message if failed

            if nargin < 4
                n_samples = 50;
            end
            if nargin < 5
                tol = 1e-5;
            end

            % Compute reachable set
            output_sets = layer.reach(input_is, method);
            % Handle different return types: single set, array, or cell array
            if iscell(output_sets)
                % Already a cell array
            elseif length(output_sets) > 1
                % Array of ImageStars - convert to cell array
                temp = cell(1, length(output_sets));
                for idx = 1:length(output_sets)
                    temp{idx} = output_sets(idx);
                end
                output_sets = temp;
            else
                % Single ImageStar
                output_sets = {output_sets};
            end

            passed = true;
            msg = '';

            n_pred = input_is.numPred;

            if n_pred > 0
                % Test corner cases
                corners = soundness_test_utils.get_predicate_corners(...
                    n_pred, input_is.pred_lb, input_is.pred_ub);

                for i = 1:size(corners, 2)
                    alpha = corners(:, i);

                    % Check constraint satisfaction
                    if ~isempty(input_is.C) && any(input_is.C * alpha > input_is.d + tol)
                        continue;
                    end

                    % Compute concrete input
                    input_concrete = input_is.V(:,:,:,1);
                    for k = 1:n_pred
                        input_concrete = input_concrete + alpha(k) * input_is.V(:,:,:,k+1);
                    end

                    % Compute concrete output
                    output_concrete = layer.evaluate(input_concrete);

                    % Check containment
                    contained = false;
                    for j = 1:length(output_sets)
                        if soundness_test_utils.verify_imagestar_containment(output_sets{j}, output_concrete, tol)
                            contained = true;
                            break;
                        end
                    end

                    if ~contained
                        passed = false;
                        msg = sprintf('Corner %d not contained in output', i);
                        return;
                    end
                end
            end

            % Test random samples
            samples = soundness_test_utils.sample_imagestar(input_is, n_samples);

            for i = 1:n_samples
                output_concrete = layer.evaluate(samples{i});

                contained = false;
                for j = 1:length(output_sets)
                    if soundness_test_utils.verify_imagestar_containment(output_sets{j}, output_concrete, tol)
                        contained = true;
                        break;
                    end
                end

                if ~contained
                    passed = false;
                    msg = sprintf('Random sample %d not contained in output', i);
                    return;
                end
            end
        end

        function results = verify_network_soundness_layer_by_layer(nnv_net, input_set, method, n_samples, tol, verbose)
            % VERIFY_NETWORK_SOUNDNESS_LAYER_BY_LAYER Layer-by-layer soundness check
            %
            % This function propagates both:
            %   1. Concrete samples through layer.evaluate()
            %   2. Reachable sets through layer.reach()
            % and verifies that at each layer, all concrete outputs are
            % contained within the reachable set bounds.
            %
            % Input:
            %   nnv_net - NNV network (NN object with Layers cell array)
            %   input_set - ImageStar or Star input set
            %   method - reach method (default: 'approx-star')
            %   n_samples - number of random samples (default: 50)
            %   tol - numerical tolerance (default: 1e-5)
            %   verbose - print detailed output (default: true)
            %
            % Output:
            %   results - struct array with fields:
            %     .layer_idx - layer index
            %     .layer_name - layer class name
            %     .passed - true if soundness verified
            %     .num_samples_tested - number of samples tested
            %     .num_violations - number of soundness violations
            %     .eval_range - [min, max] of evaluate outputs
            %     .reach_range - [min, max] of reach bounds
            %     .violation_details - details of first violation (if any)

            if nargin < 3 || isempty(method)
                method = 'approx-star';
            end
            if nargin < 4 || isempty(n_samples)
                n_samples = 50;
            end
            if nargin < 5 || isempty(tol)
                tol = 1e-5;
            end
            if nargin < 6 || isempty(verbose)
                verbose = true;
            end

            n_layers = length(nnv_net.Layers);
            results = struct('layer_idx', {}, 'layer_name', {}, 'passed', {}, ...
                'num_samples_tested', {}, 'num_violations', {}, ...
                'eval_range', {}, 'reach_range', {}, 'violation_details', {});

            if verbose
                fprintf('\n=== Layer-by-Layer Soundness Verification ===\n');
                fprintf('Network: %d layers\n', n_layers);
                fprintf('Method: %s\n', method);
                fprintf('Samples: %d\n\n', n_samples);
            end

            % Generate sample inputs from input set
            is_imagestar = isa(input_set, 'ImageStar');
            if is_imagestar
                samples = soundness_test_utils.sample_imagestar(input_set, n_samples);
            else
                samples_mat = soundness_test_utils.sample_star(input_set, n_samples);
                samples = cell(1, n_samples);
                for i = 1:n_samples
                    samples{i} = samples_mat(:, i);
                end
            end

            % Initialize current inputs
            current_samples = samples;
            current_reach = input_set;

            for layer_idx = 1:n_layers
                layer = nnv_net.Layers{layer_idx};
                layer_name = class(layer);

                if verbose
                    fprintf('Layer %d: %s... ', layer_idx, layer_name);
                end

                result = struct();
                result.layer_idx = layer_idx;
                result.layer_name = layer_name;
                result.num_samples_tested = n_samples;
                result.num_violations = 0;
                result.violation_details = {};

                try
                    % Propagate samples through evaluate
                    new_samples = cell(1, n_samples);
                    eval_outputs = [];
                    for i = 1:n_samples
                        new_samples{i} = layer.evaluate(current_samples{i});
                        eval_outputs = [eval_outputs; new_samples{i}(:)];
                    end
                    result.eval_range = [min(eval_outputs), max(eval_outputs)];

                    % Propagate reach set
                    if isa(current_reach, 'ImageStar')
                        new_reach = layer.reach(current_reach, method);
                    elseif isa(current_reach, 'Star')
                        new_reach = layer.reach(current_reach, method);
                    else
                        % Cell array of sets - take union
                        new_reach_cell = {};
                        for r = 1:length(current_reach)
                            R_r = layer.reach(current_reach{r}, method);
                            if iscell(R_r)
                                new_reach_cell = [new_reach_cell, R_r];
                            elseif length(R_r) > 1
                                for rr = 1:length(R_r)
                                    new_reach_cell{end+1} = R_r(rr);
                                end
                            else
                                new_reach_cell{end+1} = R_r;
                            end
                        end
                        new_reach = new_reach_cell;
                    end

                    % Normalize reach set to cell array
                    if iscell(new_reach)
                        reach_sets = new_reach;
                    elseif length(new_reach) > 1
                        reach_sets = cell(1, length(new_reach));
                        for idx = 1:length(new_reach)
                            reach_sets{idx} = new_reach(idx);
                        end
                    else
                        reach_sets = {new_reach};
                    end

                    % Get reach bounds
                    reach_lb = inf;
                    reach_ub = -inf;
                    for r = 1:length(reach_sets)
                        R = reach_sets{r};
                        if isa(R, 'ImageStar')
                            [lb, ub] = R.getRanges();
                            reach_lb = min(reach_lb, min(lb(:)));
                            reach_ub = max(reach_ub, max(ub(:)));
                        elseif isa(R, 'Star')
                            n_dim = R.dim;
                            for d = 1:n_dim
                                lb_d = R.getMin(d, 'linprog');
                                ub_d = R.getMax(d, 'linprog');
                                reach_lb = min(reach_lb, lb_d);
                                reach_ub = max(reach_ub, ub_d);
                            end
                        end
                    end
                    result.reach_range = [reach_lb, reach_ub];

                    % Verify containment of each sample using bounds check
                    % (simpler than full set containment, appropriate for approx methods)
                    for i = 1:n_samples
                        output_concrete = new_samples{i};

                        % Get elementwise bounds for comparison
                        [elem_lb, elem_ub] = soundness_test_utils.get_elementwise_bounds(reach_sets);

                        % Check if concrete output is within bounds
                        out_vec = output_concrete(:);
                        if length(out_vec) == length(elem_lb)
                            in_bounds = all(out_vec >= elem_lb - tol) && all(out_vec <= elem_ub + tol);
                        else
                            % Shape mismatch - use scalar bounds
                            in_bounds = min(out_vec) >= reach_lb - tol && max(out_vec) <= reach_ub + tol;
                        end

                        if ~in_bounds
                            result.num_violations = result.num_violations + 1;
                            if isempty(result.violation_details)
                                % Store first violation details
                                result.violation_details = struct(...
                                    'sample_idx', i, ...
                                    'output_range', [min(output_concrete(:)), max(output_concrete(:))], ...
                                    'reach_range', result.reach_range);
                            end
                        end
                    end

                    result.passed = (result.num_violations == 0);

                    if verbose
                        if result.passed
                            fprintf('PASSED (eval: [%.4f, %.4f], reach: [%.4f, %.4f])\n', ...
                                result.eval_range(1), result.eval_range(2), ...
                                result.reach_range(1), result.reach_range(2));
                        else
                            fprintf('FAILED (%d violations)\n', result.num_violations);
                            fprintf('         Eval range:  [%.6f, %.6f]\n', result.eval_range(1), result.eval_range(2));
                            fprintf('         Reach range: [%.6f, %.6f]\n', result.reach_range(1), result.reach_range(2));
                        end
                    end

                    % Update for next layer
                    current_samples = new_samples;
                    current_reach = new_reach;

                catch ME
                    result.passed = false;
                    result.num_violations = -1;
                    result.violation_details = struct('error', ME.message);
                    result.eval_range = [NaN, NaN];
                    result.reach_range = [NaN, NaN];

                    if verbose
                        fprintf('ERROR: %s\n', ME.message);
                    end
                end

                results(layer_idx) = result;
            end

            if verbose
                fprintf('\n=== Summary ===\n');
                num_passed = sum([results.passed]);
                num_failed = sum(~[results.passed]);
                fprintf('Passed: %d/%d layers\n', num_passed, n_layers);
                if num_failed > 0
                    fprintf('Failed layers:\n');
                    for i = 1:n_layers
                        if ~results(i).passed
                            fprintf('  Layer %d (%s): %d violations\n', ...
                                results(i).layer_idx, results(i).layer_name, ...
                                results(i).num_violations);
                        end
                    end
                end
            end
        end

        function [lb, ub] = get_output_bounds(output_set, tol)
            % GET_OUTPUT_BOUNDS Extract bounds from ImageStar or Star output
            %
            % Input:
            %   output_set - ImageStar, Star, or cell array of sets
            %   tol - tolerance (unused, for API consistency)
            %
            % Output:
            %   lb, ub - lower/upper bound vectors

            if nargin < 2
                tol = 1e-6;
            end

            % Normalize to cell array
            if iscell(output_set)
                sets = output_set;
            elseif length(output_set) > 1
                sets = cell(1, length(output_set));
                for i = 1:length(output_set)
                    sets{i} = output_set(i);
                end
            else
                sets = {output_set};
            end

            % Initialize bounds
            lb = [];
            ub = [];

            for i = 1:length(sets)
                S = sets{i};
                if isa(S, 'ImageStar')
                    [lb_i, ub_i] = S.getRanges();
                    lb_i = lb_i(:);
                    ub_i = ub_i(:);
                elseif isa(S, 'Star')
                    n = S.dim;
                    lb_i = zeros(n, 1);
                    ub_i = zeros(n, 1);
                    for j = 1:n
                        lb_i(j) = S.getMin(j, 'linprog');
                        ub_i(j) = S.getMax(j, 'linprog');
                    end
                else
                    continue;
                end

                if isempty(lb)
                    lb = lb_i;
                    ub = ub_i;
                else
                    lb = min(lb, lb_i);
                    ub = max(ub, ub_i);
                end
            end
        end

        function in_bounds = check_point_in_bounds(point, lb, ub, tol)
            % CHECK_POINT_IN_BOUNDS Check if point is within bounds
            %
            % Input:
            %   point - vector to check
            %   lb, ub - lower/upper bounds
            %   tol - tolerance (default: 1e-6)
            %
            % Output:
            %   in_bounds - true if point is within [lb-tol, ub+tol]

            if nargin < 4
                tol = 1e-6;
            end

            point = point(:);
            lb = lb(:);
            ub = ub(:);

            in_bounds = all(point >= lb - tol) && all(point <= ub + tol);
        end

        function [lb, ub] = get_elementwise_bounds(reach_sets)
            % GET_ELEMENTWISE_BOUNDS Get elementwise bounds from reach sets
            %
            % Input:
            %   reach_sets - cell array of ImageStar or Star sets
            %
            % Output:
            %   lb, ub - elementwise lower/upper bounds as vectors

            lb = [];
            ub = [];

            for i = 1:length(reach_sets)
                S = reach_sets{i};
                if isa(S, 'ImageStar')
                    [lb_i, ub_i] = S.getRanges();
                    lb_i = lb_i(:);
                    ub_i = ub_i(:);
                elseif isa(S, 'Star')
                    n = S.dim;
                    lb_i = zeros(n, 1);
                    ub_i = zeros(n, 1);
                    for j = 1:n
                        lb_i(j) = S.getMin(j, 'linprog');
                        ub_i(j) = S.getMax(j, 'linprog');
                    end
                else
                    continue;
                end

                if isempty(lb)
                    lb = lb_i;
                    ub = ub_i;
                else
                    % Take union of bounds
                    if length(lb_i) == length(lb)
                        lb = min(lb, lb_i);
                        ub = max(ub, ub_i);
                    end
                end
            end
        end

    end
end
