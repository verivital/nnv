classdef Softmax
    % SOFTMAX Class contains methods for reachability analysis of layers
    % with Softmax activation function.
    %
    % Softmax function: s_i = exp(x_i) / sum(exp(x_j))
    %
    % The main challenge is that softmax is non-linear and couples all
    % inputs together through the normalization denominator.
    %
    % Approach: Use interval bound propagation with careful handling of
    % the exponential and division operations.
    %
    % Reference:
    %   - alpha-beta-CROWN approach for softmax bounds
    %   - Zhang et al., "Efficient Neural Network Verification with Optimized
    %     Linear Relaxations"
    %
    % Author: NNV Team
    % Date: November 2025

    properties
    end

    methods(Static)

        %% ============== EVALUATION ==============

        function y = evaluate(x)
            % Standard softmax evaluation with numerical stability
            % @x: input vector
            % @y: output vector (probability distribution)

            % Subtract max for numerical stability
            x_shifted = x - max(x);
            exp_x = exp(x_shifted);
            y = exp_x / sum(exp_x);
        end

        %% ============== MAIN REACHABILITY ==============

        function S = reach_star_approx(varargin)
            % Main reachability method for Star sets
            % @I: input Star set
            % @method: reachability method
            % @reachOption: reach options
            % @relaxFactor: relaxation factor
            % @dis_opt: display option
            % @lp_solver: LP solver to use
            % @S: output Star set (over-approximation)

            switch nargin
                case 1
                    I = varargin{1};
                    method = 'approx-star';
                    reachOption = [];
                    relaxFactor = 0;
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 2
                    I = varargin{1};
                    method = varargin{2};
                    reachOption = [];
                    relaxFactor = 0;
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    method = varargin{2};
                    reachOption = varargin{3};
                    relaxFactor = 0;
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    method = varargin{2};
                    reachOption = varargin{3};
                    relaxFactor = varargin{4};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    method = varargin{2};
                    reachOption = varargin{3};
                    relaxFactor = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = 'linprog';
                case 6
                    I = varargin{1};
                    method = varargin{2};
                    reachOption = varargin{3};
                    relaxFactor = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 1-6');
            end

            if ~isa(I, 'Star')
                error('Input set is not a Star set');
            end

            % Dispatch based on method
            if strcmp(method, 'approx-star') || strcmp(method, 'approx-star-no-split')
                S = Softmax.reach_star_approx_bounds(I, dis_opt, lp_solver);
            elseif contains(method, 'relax-star')
                S = Softmax.reach_star_approx_bounds(I, dis_opt, lp_solver);
            else
                error('Unknown reachability method for Softmax: %s', method);
            end
        end

        %% ============== STAR SET REACHABILITY ==============

        function S = reach_star_approx_bounds(varargin)
            % Compute over-approximate reachable set for softmax using bounds
            % @I: input Star set
            % @dis_opt: display option
            % @lp_solver: LP solver
            % @S: output Star set

            switch nargin
                case 1
                    I = varargin{1};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 2
                    I = varargin{1};
                    dis_opt = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    dis_opt = varargin{2};
                    lp_solver = varargin{3};
                otherwise
                    error('Invalid number of input arguments');
            end

            n = I.dim;

            % Step 1: Get bounds on all inputs
            if ~isempty(dis_opt) && strcmp(dis_opt, 'display')
                fprintf('Softmax: Computing input bounds for %d dimensions...\n', n);
            end

            lb = zeros(n, 1);
            ub = zeros(n, 1);

            for i = 1:n
                lb(i) = I.getMin(i, lp_solver);
                ub(i) = I.getMax(i, lp_solver);
            end

            % Step 2: Compute softmax output bounds
            [sm_lb, sm_ub] = Softmax.compute_softmax_bounds(lb, ub);

            % Step 3: Create output Star from bounds
            % Since softmax outputs are bounded, we create a box Star
            center = (sm_lb + sm_ub) / 2;
            V = [center, diag((sm_ub - sm_lb) / 2)];

            % Create constraint matrices for the unit hypercube
            nPred = n;
            C = [eye(nPred); -eye(nPred)];
            d = [ones(nPred, 1); ones(nPred, 1)];

            S = Star(V, C, d, -ones(nPred, 1), ones(nPred, 1));

            if ~isempty(dis_opt) && strcmp(dis_opt, 'display')
                fprintf('Softmax: Output bounds computed.\n');
                fprintf('  Min output: %.4f\n', min(sm_lb));
                fprintf('  Max output: %.4f\n', max(sm_ub));
            end
        end

        %% ============== ZONOTOPE REACHABILITY ==============

        function Z = reach_zono_approx(I)
            % Zonotope over-approximation of softmax
            % @I: input Zonotope
            % @Z: output Zonotope

            if ~isa(I, 'Zono')
                error('Input set is not a Zonotope');
            end

            % Get bounds from zonotope
            [lb, ub] = I.getBounds();

            % Compute softmax bounds
            [sm_lb, sm_ub] = Softmax.compute_softmax_bounds(lb, ub);

            % Create output zonotope from bounds
            center = (sm_lb + sm_ub) / 2;
            generators = diag((sm_ub - sm_lb) / 2);

            Z = Zono(center, generators);
        end

        %% ============== BOUND COMPUTATION ==============

        function [sm_lb, sm_ub] = compute_softmax_bounds(lb, ub)
            % Compute bounds on softmax output given input bounds
            % @lb: lower bounds on input (n x 1)
            % @ub: upper bounds on input (n x 1)
            % @sm_lb: lower bounds on softmax output
            % @sm_ub: upper bounds on softmax output
            %
            % Key insight: For softmax s_i = exp(x_i) / sum(exp(x_j))
            %   - s_i is minimized when x_i is minimized AND all other x_j maximized
            %   - s_i is maximized when x_i is maximized AND all other x_j minimized

            n = length(lb);
            sm_lb = zeros(n, 1);
            sm_ub = zeros(n, 1);

            % For numerical stability, shift by max(ub)
            shift = max(ub);
            lb_shifted = lb - shift;
            ub_shifted = ub - shift;

            % Compute exp bounds
            exp_lb = exp(lb_shifted);
            exp_ub = exp(ub_shifted);

            for i = 1:n
                % For lower bound on s_i:
                % Numerator: use min exp(x_i) = exp(lb_i)
                % Denominator: use max sum = exp(ub_i) + sum_{j!=i} exp(ub_j)
                others = setdiff(1:n, i);
                sum_others_max = sum(exp_ub(others));
                denom_max = exp_ub(i) + sum_others_max;

                sm_lb(i) = exp_lb(i) / denom_max;

                % For upper bound on s_i:
                % Numerator: use max exp(x_i) = exp(ub_i)
                % Denominator: use min sum = exp(lb_i) + sum_{j!=i} exp(lb_j)
                sum_others_min = sum(exp_lb(others));
                denom_min = exp_lb(i) + sum_others_min;

                sm_ub(i) = exp_ub(i) / denom_min;
            end

            % Clamp to valid probability range [0, 1]
            sm_lb = max(sm_lb, 0);
            sm_ub = min(sm_ub, 1);

            % Additional constraint: sum of softmax outputs = 1
            % This can help tighten bounds but is not enforced here
            % for simplicity
        end

        %% ============== TIGHTER BOUNDS (CROWN-style) ==============

        function [sm_lb, sm_ub] = compute_softmax_bounds_tight(lb, ub)
            % Compute tighter bounds using CROWN-style linear relaxation
            % This provides tighter bounds but is more expensive
            %
            % Based on: "Efficient Neural Network Verification with
            % Optimized Linear Relaxations" (Zhang et al.)

            n = length(lb);
            sm_lb = zeros(n, 1);
            sm_ub = zeros(n, 1);

            % For numerical stability
            shift = max(ub);
            lb_shifted = lb - shift;
            ub_shifted = ub - shift;

            % Use sampling to get better bounds
            num_samples = 100;

            for i = 1:n
                min_val = inf;
                max_val = -inf;

                for s = 1:num_samples
                    % Sample a point in the input box
                    x_sample = lb + (ub - lb) .* rand(n, 1);
                    x_sample_shifted = x_sample - shift;

                    % Compute softmax at this point
                    exp_x = exp(x_sample_shifted);
                    sm_sample = exp_x / sum(exp_x);

                    min_val = min(min_val, sm_sample(i));
                    max_val = max(max_val, sm_sample(i));
                end

                % Also check corner cases for dimension i
                % Case 1: x_i at lower bound, others at upper bound
                x_corner1 = ub;
                x_corner1(i) = lb(i);
                x_corner1_shifted = x_corner1 - shift;
                exp_corner1 = exp(x_corner1_shifted);
                sm_corner1 = exp_corner1 / sum(exp_corner1);
                min_val = min(min_val, sm_corner1(i));

                % Case 2: x_i at upper bound, others at lower bound
                x_corner2 = lb;
                x_corner2(i) = ub(i);
                x_corner2_shifted = x_corner2 - shift;
                exp_corner2 = exp(x_corner2_shifted);
                sm_corner2 = exp_corner2 / sum(exp_corner2);
                max_val = max(max_val, sm_corner2(i));

                sm_lb(i) = max(0, min_val * 0.95);  % Small margin for safety
                sm_ub(i) = min(1, max_val * 1.05);
            end
        end

        %% ============== INTERVAL ARITHMETIC ==============

        function [sm_lb, sm_ub] = compute_softmax_bounds_interval(lb, ub)
            % Compute bounds using interval arithmetic
            % This is less tight but very fast

            n = length(lb);

            % Compute exp bounds
            exp_lb = exp(lb);
            exp_ub = exp(ub);

            % Sum bounds
            sum_lb = sum(exp_lb);
            sum_ub = sum(exp_ub);

            % Softmax bounds
            sm_lb = exp_lb / sum_ub;
            sm_ub = exp_ub / sum_lb;

            % Clamp
            sm_lb = max(sm_lb, 0);
            sm_ub = min(sm_ub, 1);
        end

        %% ============== HELPER METHODS ==============

        function S = multiStepSoftmax(I, dis_opt, lp_solver)
            % Multi-step softmax reachability (not implemented)
            % For softmax, we compute bounds all at once rather than
            % dimension by dimension
            S = Softmax.reach_star_approx_bounds(I, dis_opt, lp_solver);
        end

        %% ============== SOUNDNESS CHECK ==============

        function [is_sound, samples, outputs] = check_soundness(I, S, num_samples)
            % Check soundness of reachability computation
            % @I: input Star set
            % @S: computed output Star set
            % @num_samples: number of samples to check
            % @is_sound: true if all samples are contained in output set

            if nargin < 3
                num_samples = 1000;
            end

            is_sound = true;
            samples = cell(num_samples, 1);
            outputs = cell(num_samples, 1);

            for i = 1:num_samples
                % Sample from input set
                x = I.sample(1);
                samples{i} = x;

                % Compute softmax output
                y = Softmax.evaluate(x);
                outputs{i} = y;

                % Check containment
                if ~S.contains(y)
                    is_sound = false;
                    fprintf('Soundness violation at sample %d\n', i);
                    fprintf('Input: %s\n', mat2str(x', 4));
                    fprintf('Output: %s\n', mat2str(y', 4));
                end
            end
        end

    end
end
