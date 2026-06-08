classdef ElementwiseAffineLayer < handle
    % The ElementwiseAffineLayer layer class in CNN
    % author: Neelanjana Pal
    % date: 6/28/2021
    % update: Diego Manzanas
    %              10/31/2022
    %              Refactor properies and evaluation to mirror nnet.onnx.layer.ElementwiseAffineLayer
    
    properties
        Name = 'elementwise_affine_layer';
        Scale
        Offset
        DoScale
        DoOffset
    end
    
    
    % constructor
    methods
        
        % constructor of the class
        function obj = ElementwiseAffineLayer(varargin)    
            
            switch nargin
                
                case 5
                    obj.Name = varargin{1};
                    obj.Scale = varargin{2};
                    obj.Offset = varargin{3};
                    obj.DoScale = varargin{4};
                    obj.DoOffset = varargin{5};
                case 4
                    obj.Scale = varargin{1};
                    obj.Offset = varargin{2};
                    obj.DoScale = varargin{3};
                    obj.DoOffset = varargin{4};
                otherwise
                    error('Invalid number of inputs (should be 4 or 5)');
            end 
             
        end
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(obj, x)
            % evaluate elementwise affine layer
            y = x;
            if obj.DoScale
                s = obj.Scale;
                if ~isscalar(s), s = align_to_input(s, size(x)); end
                y = y .* s;
            end
            if obj.DoOffset
                o = obj.Offset;
                if ~isscalar(o), o = align_to_input(o, size(x)); end
                y = y + o;
            end
        end
       
    end   
     
    methods % reachability method
        
        %(reachability analysis using imagestar)
        function image = reach_star_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, "Star")
                error('Input set is not an ImageStar or Star');
            end
            
            if isa(in_image, "ImageStar")
                V = in_image.V; % Initialize value dimensions
    %             n = in_image.numPred;
                % Process scaling first
                if obj.DoScale
                    if isscalar(obj.Scale)
                        V = double(obj.Scale)*in_image.V;
                    elseif length(size(obj.Scale)) == length(size(in_image.V)) && size(obj.Scale) == size(in_image.V)
                        V = double(obj.Scale).*in_image;
                    elseif length(obj.Scale) == size(in_image.V,1)
                        for k=1:length(obj.Scale)
                            V(k, : , : , :) = V(k, : , : , :) * obj.Scale(k);
                        end
                    elseif length(obj.Scale) == size(in_image.V,2)
                        for k=1:length(obj.Scale)
                            V(:, k , : , :) = V(:, k , : , :) * obj.Scale(k);
                        end
                    elseif length(obj.Scale) == size(in_image.V,3)
                        for k=1:length(obj.Scale)
                            V(:, : , k , :) = V(:, : , k , :) * obj.Scale(k);
                        end
                    elseif length(obj.Scale) == size(in_image.V,4)
                        for k=1:length(obj.Scale)
                            V(:, : , : , k) = V(:, : , : , k) * obj.Scale(k);
                        end
                    else
                        % Fallback: try broadcast-shape align like evaluate.
                        % numel-aware reshape that lands the non-singleton dims
                        % of Scale onto matching V spatial dims.
                        sz_v = size(V);
                        sz_v(end) = 1;   % treat the last (predicate) dim as 1
                        try
                            s_aligned = align_to_input(obj.Scale, sz_v);
                            V = V .* s_aligned;
                        catch
                            error('TODO: add support for other tensor shapes (Scale)')
                        end
                    end
                end
    
                % Process offset last
                if obj.DoOffset
                    % offset = squeeze(obj.Offset);
    %                 a = size(offset);
                    if isscalar(obj.Offset)
                        V = V + obj.Offset;
                    elseif length(obj.Offset) == size(in_image.V,1)
                        for k=1:length(obj.Offset)
                            V(k, : , : , :) = V(k, : , : , :) + obj.Offset(k);
                        end
                    elseif length(obj.Offset) == size(in_image.V,2)
                        for k=1:length(obj.Offset)
                            V(:, k , : , :) = V(:, k , : , :) + obj.Offset(k);
                        end
                    elseif length(obj.Offset) == size(in_image.V,3)
                        for k=1:length(obj.Offset)
                            V(:, : , k , 1) = V(:, : , k , 1) + obj.Offset(k); % bias only added to first column
                        end
                    elseif length(obj.Offset) == size(in_image.V,4)
                        for k=1:length(obj.Offset)
                            V(:, : , : , k) = V(:, : , : , k) + obj.Offset(k);
                        end
                    else
                        % Fallback: align bias shape to V's spatial dims via
                        % align_to_input. ONLY add the bias to V(:,:,:,1) (the
                        % center column) — basis vectors should not get
                        % shifted, otherwise predicate semantics are wrong.
                        sz_v_spatial = size(V); sz_v_spatial = sz_v_spatial(1:max(1, ndims(V)-1));
                        try
                            o_aligned = align_to_input(obj.Offset, sz_v_spatial);
                            V(:,:,:,1) = V(:,:,:,1) + o_aligned;
                        catch
                            error('TODO: add support for other tensor shapes (Offset).')
                        end
                    end
                end
    
                % return output set
                image = ImageStar(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            else % reachability with Star sets
                % In general, this layer is used as the bias addition following a FullyConnectedLayer
                % Scale (weights)
                image = in_image; % create copy to perform operations (if need to)
                if obj.DoScale
                    s = obj.Scale;
                    if isscalar(s)
                        % Broadcast scalar scale over all dimensions:
                        % use a diagonal matrix of size [dim, dim]
                        image = image.affineMap(diag(double(s) * ones(image.dim, 1)), []);
                    elseif numel(s) == image.dim
                        % ND scale collapses to a per-dim diagonal scaling
                        image = image.affineMap(diag(double(s(:))), []);
                    else
                        image = image.affineMap(s, []); % full matrix
                    end
                end
                % Offset (bias)
                if obj.DoOffset
                    o = obj.Offset;
                    if isscalar(o)
                        o = double(o) * ones(image.dim, 1);
                    elseif numel(o) == image.dim
                        o = double(o(:));   % flatten ND bias to [dim,1]
                    elseif numel(o) > image.dim && mod(numel(o), image.dim) == 0
                        % Replicated bias: e.g. tile(thresholds, k) where the
                        % flatten produced k*dim elements but the layer only
                        % takes the first dim slice. Take first dim entries.
                        flat = double(o(:));
                        o = flat(1:image.dim);
                    elseif numel(o) < image.dim && mod(image.dim, numel(o)) == 0
                        % Bias is a per-channel constant that should tile to
                        % fill image.dim. Repeat to match.
                        flat = double(o(:));
                        o = repmat(flat, image.dim / numel(flat), 1);
                    else
                        % Last-resort: pad with zeros or truncate
                        flat = double(o(:));
                        if numel(flat) < image.dim
                            flat = [flat; zeros(image.dim - numel(flat), 1)];
                        end
                        o = flat(1:image.dim);
                    end
                    image = image.affineMap(diag(ones(1,image.dim)), o); % x + b
                end
            end
            
        end
        
        % handle multiple inputs
        function S = reach_star_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
            n = length(inputs);

            % Initialize output variables
            if isa(inputs, "ImageStar")
                S(n) = ImageStar;
            else
                S(n) = Star;
            end

            % Begin computing reachability one set at a time
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_star_single_input(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_star_single_input(inputs(i));
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
            
        end
        
        
        % reachability analysis with multiple inputs
        function IS = reach(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
            
            % author: Neelanjana Pal
            % date: 6/28/2021
           
             
            switch nargin
                
                 case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                    % lp_solver = varargin{7}; do not use
                
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                
                case 4
                    obj = varargin{1};
                    in_images = varargin{2}; 
                    method = varargin{3};
                    option = varargin{4}; % computation option

                case 3
                    obj = varargin{1};
                    in_images = varargin{2}; % don't care the rest inputs
                    method = varargin{3};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5 or 6)');
            end
            
            if obj.DoScale || obj.DoOffset
                if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || contains(method, "relax-star")
                    IS = obj.reach_star_multipleInputs(in_images, option);
                else
                    error('Unsupported/Unknown reachability method, only star based methods are supported');
                end
            else
                IS = in_images;
            end
  
        end
        
    end

    methods % helper functions

        % change params to gpuArrays
        function obj = toGPU(obj)
            obj.Scale = gpuArray(obj.Scale);
            obj.Offset = gpuArray(obj.Offset);
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, "double")
                obj.Offset = double(obj.Offset);
                obj.Scale = double(obj.Scale);
            elseif strcmp(precision, "single")
                obj.Offset = single(obj.Offset);
                obj.Scale = single(obj.Scale);
            else
                error("Parameter numerical precision must be 'single' or 'double'");
            end
        end

    end
    
    
    methods(Static)
        
        % parse a trained elementwise affine layer from matlab
        function L = parse(elementwise_affine_layer)
            % @elementwise_affine_layer: a elementwise affine layer from matlab deep
            % neural network tool box
            % @L : a ElementwiseAffineLayer obj for reachability analysis purpose            
            
            if ~isa(elementwise_affine_layer, 'nnet.onnx.layer.ElementwiseAffineLayer')
                error('Input is not a Matlab nnet.onnx.layer.ElementwiseAffineLayer class');
            end
                        
            L = ElementwiseAffineLayer(elementwise_affine_layer.Name, elementwise_affine_layer.Scale, elementwise_affine_layer.Offset, ...
                elementwise_affine_layer.DoScale, elementwise_affine_layer.DoOffset);

        end

    end


end

function out = padshape_left(arr, target_ndims)
% Align arr's trailing dims with the target. Strip leading singleton dims
% if too many; prepend if too few. So bias [1,1,1,5] aligned to ndims=2
% becomes [1,5]; bias [5,48] aligned to ndims=3 becomes [1,5,48].
sz = size(arr);
% Strip leading singletons (only if they're truly singleton)
while numel(sz) > target_ndims && sz(1) == 1
    sz = sz(2:end);
end
% Prepend if still too few dims
need = target_ndims - numel(sz);
if need > 0
    sz = [ones(1, need), sz];
end
out = reshape(arr, sz);
end

function out = align_to_input(arr, x_shape)
% Smart reshape of an ONNX broadcasting bias/scale to match a given input
% shape. Squeezes arr to its non-singleton sizes, then maps each non-singleton
% bias dim to the matching input dim (by size). Falls back to padshape_left
% if no clean match is possible.
%
% Examples:
%   arr=[1,1,1,5], x_shape=[5,1] -> [5,1]
%   arr=[5,48],    x_shape=[1,5,48] -> [1,5,48]
%   arr=[1,1,16],  x_shape=[28,28,16] -> [1,1,16]
sz_a = size(arr);
ns_a = sz_a(sz_a > 1);   % non-singleton dims of arr (in order)
if isempty(ns_a)
    out = arr; return;   % all-singleton, broadcast-compatible
end

% If arr has exactly one non-singleton dim, place it where input has the same size
if numel(ns_a) == 1
    target_size = ones(1, max(2, numel(x_shape)));
    placed = false;
    val = ns_a(1);
    for k = 1:numel(x_shape)
        if x_shape(k) == val && ~placed
            target_size(k) = val; placed = true;
        end
    end
    if placed
        out = reshape(arr, target_size); return;
    end
end

% If non-singleton sizes are a subsequence of x_shape, place each at matching position
target_size = ones(1, max(2, numel(x_shape)));
remaining = ns_a;
ok = true;
for k = 1:numel(x_shape)
    if ~isempty(remaining) && x_shape(k) == remaining(1)
        target_size(k) = remaining(1);
        remaining = remaining(2:end);
    end
end
if isempty(remaining) && ok
    out = reshape(arr, target_size); return;
end

% Last resort: keep arr's shape, padshape_left handles dim mismatch
out = padshape_left(arr, numel(x_shape));
end

