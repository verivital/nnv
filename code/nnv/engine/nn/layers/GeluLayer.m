classdef GeluLayer < handle
    % GeluLayer - GELU (Gaussian Error Linear Unit) activation layer
    %
    % GELU(x) = x * Φ(x) where Φ is the CDF of standard normal distribution
    % Approximation: GELU(x) ≈ 0.5 * x * (1 + tanh(√(2/π) * (x + 0.044715 * x³)))
    %
    % GELU is commonly used in Transformer architectures (BERT, GPT, ViT).
    %
    % Author: NNV Team
    % Date: December 2025

    properties
        Name = 'gelu_layer';
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
        % Which GELU formula this layer models. MATLAB's geluLayer / the ONNX
        % Gelu op default to the EXACT erf form ('none'); NNV historically used
        % the tanh approximation. evaluate() AND the reach bounds must use the
        % SAME formula as the source network, else the bounds (computed for one
        % form) need not contain the true outputs of the other -- the two forms
        % differ by up to ~4.7e-4 [2][25]. parse() syncs this to the source.
        Approximation = 'tanh';   % 'tanh' (approx) or 'none' (exact erf-GELU)
    end

    methods

        %% Constructor
        function obj = GeluLayer(varargin)
            % Constructor for GeluLayer
            %
            % Usage:
            %   layer = GeluLayer()
            %   layer = GeluLayer(name)
            %   layer = GeluLayer(name, numInputs, inputNames, numOutputs, outputNames)

            switch nargin
                case 0
                    % Default
                case 1
                    obj.Name = varargin{1};
                case 5
                    obj.Name = varargin{1};
                    obj.NumInputs = varargin{2};
                    obj.InputNames = varargin{3};
                    obj.NumOutputs = varargin{4};
                    obj.OutputNames = varargin{5};
                otherwise
                    error('Invalid number of arguments');
            end
        end

        %% Evaluation
        function y = evaluate(obj, input)
            % Evaluate GELU activation using the layer's configured formula
            % (exact erf-GELU when Approximation='none', tanh approximation
            % otherwise). gelu_apply is the single source of truth shared with
            % the reach bounds so reach is sound w.r.t. evaluate.
            y = obj.gelu_apply(input);
        end

        function y = gelu_apply(obj, x)
            % The actual GELU formula, selected by obj.Approximation. Used by
            % both evaluate() and the reach bound endpoints (gelu_scalar) so the
            % two can never diverge.
            if strcmpi(obj.Approximation, 'none')
                y = 0.5 .* x .* (1 + erf(x ./ sqrt(2)));            % exact GELU
            else
                s2p = sqrt(2/pi);
                y = 0.5 .* x .* (1 + tanh(s2p .* (x + 0.044715 .* x.^3)));  % tanh approx
            end
        end

        function y = evaluateSequence(obj, input)
            % Evaluate for sequence data
            y = obj.evaluate(input);
        end

        %% Reachability Analysis
        function images = reach_star_single_input(obj, in_image, method, relaxFactor, dis_opt, lp_solver)
            % Reachability for single ImageStar or Star input
            %
            % Uses bounds-based approximation for GELU.

            if nargin < 6, lp_solver = 'linprog'; end
            if nargin < 5, dis_opt = []; end
            if nargin < 4, relaxFactor = 0; end
            if nargin < 3, method = 'approx-star'; end

            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star')
                error('Input must be ImageStar or Star');
            end

            if isa(in_image, 'ImageStar')
                % Get input bounds
                if ~isempty(in_image.im_lb) && ~isempty(in_image.im_ub)
                    lb = in_image.im_lb;
                    ub = in_image.im_ub;
                else
                    [lb, ub] = in_image.estimateRanges();
                end

                % Compute GELU bounds element-wise
                [out_lb, out_ub] = obj.compute_gelu_bounds(lb, ub);

                % Create output ImageStar
                images = ImageStar(out_lb, out_ub);

            else % Star
                n = in_image.dim;
                lb = zeros(n, 1);
                ub = zeros(n, 1);

                for i = 1:n
                    lb(i) = in_image.getMin(i, lp_solver);
                    ub(i) = in_image.getMax(i, lp_solver);
                end

                % Compute GELU bounds
                [out_lb, out_ub] = obj.compute_gelu_bounds(lb, ub);

                % Create output Star from bounds
                images = Star(out_lb, out_ub);
            end
        end

        function image = reach_zono(obj, in_image)
            % Reachability for ImageZono or Zono

            if ~isa(in_image, 'ImageZono') && ~isa(in_image, 'Zono')
                error('Input must be ImageZono or Zono');
            end

            if isa(in_image, 'ImageZono')
                % ImageZono uses getRanges(), not getBounds()
                [lb, ub] = in_image.getRanges();
                [out_lb, out_ub] = obj.compute_gelu_bounds(lb, ub);
                image = ImageZono(out_lb, out_ub);
            else
                % Zono has getBounds()
                [lb, ub] = in_image.getBounds();
                [out_lb, out_ub] = obj.compute_gelu_bounds(lb, ub);
                center = (out_lb + out_ub) / 2;
                generators = diag((out_ub - out_lb) / 2);
                image = Zono(center, generators);
            end
        end

        function [out_lb, out_ub] = compute_gelu_bounds(obj, lb, ub)
            % Compute output bounds for GELU given input bounds
            %
            % GELU is approximately monotonic for x > -1.5
            % For x < -1.5, GELU(x) ≈ 0

            out_lb = zeros(size(lb));
            out_ub = zeros(size(ub));

            % Sound bracket around the single GELU minimum. The argmin differs
            % by variant -- erf ~ -0.75179, tanh ~ -0.75246 -- so a single split
            % point would mis-classify the monotone branch and miss the true min
            % for the other variant. Outside [x_min_lo, x_min_hi] BOTH variants
            % are strictly monotone (verified numerically): decreasing to the
            % left, increasing to the right. gelu_min is a conservative lower
            % bound below both true minima (-0.16997 / -0.17004).
            x_min_lo = -0.78;
            x_min_hi = -0.72;
            gelu_min = -0.1701;

            for i = 1:numel(lb)
                l = lb(i);
                u = ub(i);

                % GELU values at bounds (using THIS layer's variant)
                gelu_l = obj.gelu_scalar(l);
                gelu_u = obj.gelu_scalar(u);

                if u <= x_min_lo
                    % Entirely on the decreasing branch
                    out_lb(i) = gelu_u;
                    out_ub(i) = gelu_l;
                elseif l >= x_min_hi
                    % Entirely on the increasing branch
                    out_lb(i) = gelu_l;
                    out_ub(i) = gelu_u;
                else
                    % Interval may contain the minimum: min is the global min,
                    % max is at an endpoint (min is interior to the bracket)
                    out_lb(i) = gelu_min;
                    out_ub(i) = max(gelu_l, gelu_u);
                end
            end
        end

        function y = gelu_scalar(obj, x)
            % GELU at scalar/array x, using this layer's configured variant so
            % the reach bounds match evaluate().
            y = obj.gelu_apply(x);
        end

        %% Multiple inputs handling
        function images = reach_star_multipleInputs(obj, in_images, method, option, relaxFactor, dis_opt, lp_solver)
            n = length(in_images);
            images = [];

            if strcmp(option, 'parallel')
                parfor i = 1:n
                    images = [images obj.reach_star_single_input(in_images(i), method, relaxFactor, dis_opt, lp_solver)];
                end
            else
                for i = 1:n
                    images = [images obj.reach_star_single_input(in_images(i), method, relaxFactor, dis_opt, lp_solver)];
                end
            end
        end

        function images = reach_zono_multipleInputs(obj, in_images, option)
            n = length(in_images);

            if strcmp(option, 'parallel')
                parfor i = 1:n
                    images(i) = obj.reach_zono(in_images(i));
                end
            else
                for i = 1:n
                    images(i) = obj.reach_zono(in_images(i));
                end
            end
        end

        %% Main reach method
        function images = reach(varargin)
            % Main reachability method
            %
            % Usage:
            %   images = layer.reach(in_images, method, option, relaxFactor, dis_opt, lp_solver)

            switch nargin
                case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5};
                    dis_opt = varargin{6};
                    lp_solver = varargin{7};
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5};
                    dis_opt = varargin{6};
                    lp_solver = 'linprog';
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = 0;
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = 'single';
                    relaxFactor = 0;
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 2
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = 'approx-star';
                    option = 'single';
                    relaxFactor = 0;
                    dis_opt = [];
                    lp_solver = 'linprog';
                otherwise
                    error('Invalid number of arguments: %d', nargin);
            end

            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || ...
               strcmp(method, 'abs-dom') || contains(method, 'relax-star')
                images = obj.reach_star_multipleInputs(in_images, method, option, relaxFactor, dis_opt, lp_solver);
            elseif strcmp(method, 'approx-zono')
                images = obj.reach_zono_multipleInputs(in_images, option);
            else
                error('Unknown reachability method: %s', method);
            end
        end

        %% Utility methods
        function obj = toGPU(obj)
            % No parameters to move to GPU
        end

        function obj = changeParamsPrecision(obj, ~)
            % No parameters to change precision
        end

    end

    methods (Static)

        function L = parse(layer)
            % Parse MATLAB GELULayer
            %
            % Input: MATLAB nnet.cnn.layer.GELULayer
            % Output: NNV GeluLayer

            if ~isa(layer, 'nnet.cnn.layer.GELULayer')
                error('Input must be a MATLAB GELULayer');
            end

            L = GeluLayer(layer.Name, layer.NumInputs, layer.InputNames, ...
                          layer.NumOutputs, layer.OutputNames);

            % Sync the GELU formula to the source layer so evaluate/reach model
            % the ACTUAL network [2][25]. MATLAB's geluLayer exposes
            % Approximation in {"none" (exact erf, the default), "tanh"}.
            if isprop(layer, 'Approximation') && ~isempty(layer.Approximation)
                ap = lower(char(string(layer.Approximation)));
                if strcmp(ap, 'none')
                    L.Approximation = 'none';   % exact erf-GELU
                else
                    L.Approximation = 'tanh';
                end
            else
                % Property absent: MATLAB's geluLayer default is exact erf.
                L.Approximation = 'none';
            end
        end

    end

end
