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

                    % Compute the approximation error bounds
                    % Sample points to get tight bounds on the error
                    n_samples = 20;
                    x_samples = linspace(l, u, n_samples);

                    % IMPORTANT: Include the SiLU minimum point if it falls
                    % within the interval [l, u]. The SiLU minimum is a
                    % critical point that must be sampled for soundness.
                    [x_min_silu, ~] = SiLU.get_minimum();
                    if l <= x_min_silu && x_min_silu <= u
                        x_samples = [x_samples, x_min_silu];
                    end

                    error_lb = inf;
                    error_ub = -inf;

                    for j = 1:length(x_samples)
                        x_j = x_samples(j);
                        y_exact = SiLU.evaluate(x_j);
                        y_approx = y_center + dy_center * (x_j - x_center);
                        error_j = y_exact - y_approx;
                        error_lb = min(error_lb, error_j);
                        error_ub = max(error_ub, error_j);
                    end

                    % Add small tolerance for numerical robustness
                    tol = 1e-6;
                    error_lb = error_lb - tol;
                    error_ub = error_ub + tol;

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
            y_u = max(y_l_endpoint, y_u_endpoint);

            % Compute derivatives at endpoints
            dy_l = SiLU.gradient(l);
            dy_u = SiLU.gradient(u);

            % Secant line coefficients (for lower/upper bounds)
            slope = (y_u_endpoint - y_l_endpoint) / (u - l);

            % Lower bound tangent: use the smaller derivative endpoint
            if dy_l <= dy_u
                alpha_l = dy_l;
                beta_l = y_l_endpoint - dy_l * l;
            else
                alpha_l = dy_u;
                beta_l = y_u_endpoint - dy_u * u;
            end

            % Upper bound: use secant line
            alpha_u = slope;
            beta_u = y_l_endpoint - slope * l;

            % Handle the case where the function is concave in [l, u]
            % or need to adjust bounds based on convexity

            % For robustness, ensure bounds are valid
            if alpha_l > alpha_u
                % Use secant for both
                alpha_l = slope;
                beta_l = y_l_endpoint - slope * l;
            end
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
