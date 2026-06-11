classdef SiLU
    % SiLU Class contains methods for reachability analysis of layers with
    % SiLU (Sigmoid Linear Unit) activation function, also known as Swish.
    %
    % SiLU(x) = x * sigmoid(x) = x / (1 + exp(-x))
    %
    % Properties:
    %   - Smooth and differentiable everywhere
    %   - Bounded below: min(SiLU(x)) ≈ -0.278 at x ≈ -1.278
    %   - Asymptotic: SiLU(x) → x as x → +∞
    %   - Asymptotic: SiLU(x) → 0 as x → -∞
    %
    % Reference: https://arxiv.org/abs/1710.05941 (Searching for Activation Functions)
    %
    % Author: NNV Team
    % Date: November 2025

    properties
    end

    methods(Static)

        %% Evaluation Methods

        function y = evaluate(x)
            % Evaluate SiLU activation: y = x * sigmoid(x)
            % @x: input (scalar, vector, or array)
            % @y: output with same shape as input

            y = x .* (1 ./ (1 + exp(-x)));
        end

        function dy = gradient(x)
            % Compute gradient of SiLU at x
            % SiLU'(x) = sigmoid(x) + x * sigmoid(x) * (1 - sigmoid(x))
            %         = sigmoid(x) * (1 + x * (1 - sigmoid(x)))

            sig = 1 ./ (1 + exp(-x));
            dy = sig .* (1 + x .* (1 - sig));
        end

        function [x_min, y_min] = get_minimum()
            % Return the global minimum of SiLU
            % SiLU has a unique minimum at approximately x ≈ -1.278
            % where SiLU(x) ≈ -0.278

            % More precise values (computed numerically)
            x_min = -1.2784645427610737;
            y_min = -0.27846454276107373;
        end

        %% Main Reachability Method

        function S = reach(varargin)
            % Main reachability method for SiLU
            % @I: input Star set
            % @method: reachability method ('approx-star', 'approx-zono')
            % @S: output Star set(s)

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
                    error('Invalid number of input arguments');
            end

            if strcmp(method, 'approx-star') || strcmp(method, 'approx-star-no-split') || contains(method, 'relax-star')
                S = SiLU.reach_star_approx(I, dis_opt, lp_solver);
            elseif strcmp(method, 'approx-zono')
                S = SiLU.reach_zono_approx(I);
            else
                error('Unknown reachability method: %s', method);
            end
        end

        %% Star Set Reachability

        function S = reach_star_approx(varargin)
            % Approximate reachability using Star sets
            % Uses linear relaxation of SiLU

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

            if ~isa(I, 'Star')
                error('Input must be a Star set');
            end

            S = SiLU.multiStepSiLU_NoSplit(I, dis_opt, lp_solver);
        end

        function S = multiStepSiLU_NoSplit(I, dis_opt, lp_solver)
            % Multi-step SiLU reachability without splitting
            % Process all neurons in a single pass

            if ~isa(I, 'Star')
                error('Input must be a Star set');
            end

            n = I.dim;

            % Get bounds for all dimensions
            if isempty(I.predicate_lb) || isempty(I.predicate_ub)
                [lb, ub] = I.estimateRanges;
            else
                lb = I.predicate_lb;
                ub = I.predicate_ub;
            end

            % Compute bounds more precisely using LP if needed
            [lb_precise, ub_precise] = I.getRanges;

            % [19] CRASH FIX: the precise predicate-augmenting path below assumes
            % the input carries explicit predicate bounds of length nVar. A Star
            % built with the 3-arg Star(V,C,d) (halfspace intersections, some
            % NNCS/set ops) has EMPTY predicate_lb/ub, so [I.predicate_lb;
            % zeros(n,1)] would be length n (not nVar+n) and the 5-arg Star
            % constructor would error mid-network. Fall back to a SOUND
            % per-coordinate interval (box) over-approximation from the resolved
            % output ranges -- looser (drops input correlations) but correct, and
            % only taken on this uncommon input class. SiLU has a single minimum
            % (x ~ -1.278); it is decreasing left of it and increasing right.
            if isempty(I.predicate_lb) || isempty(I.predicate_ub)
                lo = lb_precise(:); hi = ub_precise(:);
                [xm, ym] = SiLU.get_minimum();
                sl = SiLU.evaluate(lo); su = SiLU.evaluate(hi);
                olb = zeros(n, 1); oub = zeros(n, 1);
                for i = 1:n
                    if hi(i) <= xm
                        olb(i) = su(i); oub(i) = sl(i);          % decreasing branch
                    elseif lo(i) >= xm
                        olb(i) = sl(i); oub(i) = su(i);          % increasing branch
                    else
                        olb(i) = ym; oub(i) = max(sl(i), su(i)); % straddles the minimum
                    end
                end
                S = Star(olb, oub);
                return;
            end

            % Initialize output
            new_V = zeros(n, I.nVar + n + 1);
            new_C = zeros(0, I.nVar + n);
            new_d = zeros(0, 1);
            new_pred_lb = [I.predicate_lb; zeros(n, 1)];
            new_pred_ub = [I.predicate_ub; zeros(n, 1)];

            % Copy existing predicate variables
            new_V(:, 1:I.nVar+1) = 0;
            new_V(:, 1) = I.V(:, 1);  % Center

            % Initialize constraint matrices from input Star
            if ~isempty(I.C)
                % Pad existing constraints with zeros for new variables
                new_C = [I.C, zeros(size(I.C, 1), n)];
                new_d = I.d;
            end

            for i = 1:n
                l = lb_precise(i);
                u = ub_precise(i);

                if l == u
                    % Constant input - exact output
                    y_val = SiLU.evaluate(l);
                    new_V(i, 1) = y_val;
                    new_pred_lb(I.nVar + i) = 0;
                    new_pred_ub(I.nVar + i) = 0;
                else
                    % Compute SiLU bounds and linear relaxation
                    [y_l, y_u, alpha_l, beta_l, alpha_u, beta_u] = SiLU.getLinearBounds(l, u);

                    % Create new predicate variable for this neuron
                    % The output is: y_i = alpha * x_i + beta + delta_i
                    % where delta_i is bounded

                    % Use the average slope for the center
                    alpha_avg = (y_u - y_l) / (u - l);
                    x_center = I.V(i, 1);
                    y_center = SiLU.evaluate(x_center);

                    % Set the center of output
                    new_V(i, 1) = y_center;

                    % Copy the linear dependency on input variables
                    % scaled by the average derivative at the center
                    dy_center = SiLU.gradient(x_center);
                    new_V(i, 2:I.nVar+1) = dy_center * I.V(i, 2:I.nVar+1);

                    % Add new predicate variable for approximation error
                    new_V(i, I.nVar + 1 + i) = 1;

                    % SOUND bound on the linearization error
                    %   e(x) = SiLU(x) - (y_center + dy_center*(x - x_center))
                    % over [l,u]. Sampling ALONE is UNSOUND: the true extrema can
                    % fall BETWEEN samples. e is C^2 with |e''| = |SiLU''| <= M
                    % (proven global max 0.5, attained at x=0), so on a subinterval
                    % of width h, e deviates from the chord of its two endpoint
                    % samples by at most M*h^2/8. Taking the per-segment chord
                    % range +/- that correction over all adjacent breakpoints is a
                    % GUARANTEED enclosure. Endpoints + the SiLU minimum + the
                    % linearization point are included as breakpoints for tightness.
                    M_silu2 = 0.5;   % proven global max |SiLU''|
                    n_samples = 20;
                    x_samples = linspace(l, u, n_samples);
                    [x_min_silu, ~] = SiLU.get_minimum();
                    if l <= x_min_silu && x_min_silu <= u
                        x_samples = [x_samples, x_min_silu];
                    end
                    if l < x_center && x_center < u
                        x_samples = [x_samples, x_center];
                    end
                    x_samples = unique(sort(x_samples(:).'));   % sorted breakpoints
                    e_samples = (x_samples .* (1 ./ (1 + exp(-x_samples)))) ...
                                - (y_center + dy_center * (x_samples - x_center));
                    error_lb = inf; error_ub = -inf;
                    for j = 1:numel(x_samples) - 1
                        h = x_samples(j+1) - x_samples(j);
                        corr = M_silu2 * h^2 / 8;                  % sound chord-deviation correction
                        error_lb = min(error_lb, min(e_samples(j), e_samples(j+1)) - corr);
                        error_ub = max(error_ub, max(e_samples(j), e_samples(j+1)) + corr);
                    end
                    error_lb = error_lb - 1e-9;   % floating-point guard only
                    error_ub = error_ub + 1e-9;

                    new_pred_lb(I.nVar + i) = error_lb;
                    new_pred_ub(I.nVar + i) = error_ub;
                end
            end

            % Create output Star
            S = Star(new_V, new_C, new_d, new_pred_lb, new_pred_ub);
        end

        function [y_l, y_u, alpha_l, beta_l, alpha_u, beta_u] = getLinearBounds(l, u)
            % Compute linear bounds for SiLU on interval [l, u]
            % Returns output bounds and linear relaxation coefficients
            %
            % Lower bound: y >= alpha_l * x + beta_l
            % Upper bound: y <= alpha_u * x + beta_u

            % Evaluate at endpoints
            y_l_endpoint = SiLU.evaluate(l);
            y_u_endpoint = SiLU.evaluate(u);

            % Get the minimum of SiLU
            [x_min, y_min] = SiLU.get_minimum();

            % Check if minimum is in the interval
            if l <= x_min && x_min <= u
                % Minimum is in interval
                y_l = y_min;
            else
                % Minimum at endpoints
                y_l = min(y_l_endpoint, y_u_endpoint);
            end

            % Upper bound is always at one of the endpoints for this range
            % (SiLU's only interior extremum on [l,u] is its minimum).
            y_u = max(y_l_endpoint, y_u_endpoint);

            % [21] SOUNDNESS: SiLU is NON-CONVEX (inflections near +/-2.4). A
            % tangent line is NOT a valid global lower bound on an interval
            % containing the dip, and the secant is an upper bound only on convex
            % sub-intervals (on concave parts it lies below the function). The old
            % code returned tangent/secant coefficients as "global" linear bounds
            % -- a silent-unsoundness trap for any caller that trusts the
            % documented contract (the current consumer, multiStepSiLU_NoSplit,
            % uses only y_l/y_u). Return SOUND -- if loose -- CONSTANT linear
            % bounds instead: y_l <= SiLU(x) <= y_u holds for all x in [l,u], so
            %   y >= 0*x + y_l   and   y <= 0*x + y_u
            % are valid global bounds. A future caller wanting a tighter sound
            % relaxation should split at the inflection/minimum first.
            alpha_l = 0; beta_l = y_l;
            alpha_u = 0; beta_u = y_u;
        end

        %% Zonotope Reachability

        function Z = reach_zono_approx(I)
            % Approximate reachability using Zonotopes

            if ~isa(I, 'Zono') && ~isa(I, 'Star')
                error('Input must be a Zono or Star');
            end

            if isa(I, 'Star')
                % Convert Star to Zono
                I = I.getZono;
            end

            c = I.c;  % Center
            V = I.V;  % Generators

            n = length(c);

            % Compute bounds
            B = I.getBox;
            lb = B.lb;
            ub = B.ub;

            % New center and generators
            new_c = zeros(n, 1);
            new_V = zeros(n, size(V, 2) + n);

            for i = 1:n
                l = lb(i);
                u = ub(i);

                if l == u
                    new_c(i) = SiLU.evaluate(l);
                    new_V(i, :) = 0;
                else
                    % Compute output range
                    y_l = SiLU.evaluate(l);
                    y_u = SiLU.evaluate(u);

                    % Check for minimum
                    [x_min, y_min] = SiLU.get_minimum();
                    if l <= x_min && x_min <= u
                        y_l_actual = y_min;
                    else
                        y_l_actual = min(y_l, y_u);
                    end
                    y_u_actual = max(y_l, y_u);

                    % Linear approximation (secant)
                    slope = (y_u - y_l) / (u - l);

                    % Apply linear transformation to existing generators
                    new_V(i, 1:size(V, 2)) = slope * V(i, :);

                    % Center of linear approximation
                    y_linear_center = y_l + slope * (c(i) - l);
                    y_exact_center = SiLU.evaluate(c(i));

                    % Compute error bound (conservative)
                    err_bound = max(abs(y_u_actual - y_exact_center), abs(y_l_actual - y_exact_center));

                    new_c(i) = y_exact_center;
                    new_V(i, size(V, 2) + i) = err_bound;  % New generator for error
                end
            end

            Z = Zono(new_c, new_V);
        end

    end
end
