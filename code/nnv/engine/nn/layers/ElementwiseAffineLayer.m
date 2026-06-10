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
            %
            % [10] FIX: a parameter whose ND SHAPE already broadcasts against x
            % (e.g. a [1,1,C] per-channel offset over an [H,W,C] feature map) must
            % use MATLAB implicit expansion -- align_to_input only looks at the
            % parameter's LENGTH and snaps its non-singleton dim onto the FIRST
            % input dim of that size, which mis-applies the [1,1,C] offset to the
            % rows when H==C. Only re-align genuinely flat/mismatched shapes
            % (e.g. a [1,1,1,N] param vs an [N,1] feature vector -- Test 21).
            y = x;
            if obj.DoScale
                s = obj.Scale;
                if ~isscalar(s) && ~elementwiseaffine_broadcastable(size(s), size(x))
                    s = align_to_input(s, size(x));
                end
                y = y .* s;
            end
            if obj.DoOffset
                o = obj.Offset;
                if ~isscalar(o) && ~elementwiseaffine_broadcastable(size(o), size(x))
                    o = align_to_input(o, size(x));
                end
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
                % Unambiguous spatial size [H,W,C]. Do NOT derive it from
                % size(V) by dropping the last dim: a ZERO-predicate (single
                % point) ImageStar has V of shape [H,W,C] (no 4th dim), so
                % size(V); spatial(end)=[] would drop the CHANNEL axis -> a
                % per-channel scale lands on the wrong dim (silent-wrong when
                % H==C, or a hard size-mismatch when H~=C).
                img_spatial = [in_image.height, in_image.width, in_image.numChannel];
                % Process scaling first.
                % [10] FIX: select the broadcast axis from the parameter's
                % SHAPE (as evaluate does), NOT its length. Length-matching the
                % non-singleton dim onto the first V dim of that size mis-places
                % a [1,1,C] channel scale onto the rows when H==C (two spatial
                % dims share a size). Scale is linear, so it multiplies the
                % center AND every basis column -- implicit expansion across the
                % trailing predicate dim handles that automatically.
                if obj.DoScale
                    if isscalar(obj.Scale)
                        V = double(obj.Scale)*in_image.V;
                    else
                        s = double(obj.Scale);
                        if ~elementwiseaffine_broadcastable(size(s), img_spatial)
                            % flat/lower-rank param (e.g. [C,1]) -> reuse the
                            % same alignment evaluate falls back to.
                            s = align_to_input(s, img_spatial);
                        end
                        V = V .* s;
                    end
                end

                % Process offset last.
                % [10] FIX: same shape-based axis selection as scale. A constant
                % bias is added ONLY to the center column V(:,:,:,1); shifting
                % the basis columns would corrupt predicate semantics and inflate
                % the set (the old dim-1/2/4 branches did exactly that, blowing
                % the reachable set up by ~|sum(a_i)| when H==C).
                if obj.DoOffset
                    if isscalar(obj.Offset)
                        V(:,:,:,1) = V(:,:,:,1) + obj.Offset;
                    else
                        o = double(obj.Offset);
                        if ~elementwiseaffine_broadcastable(size(o), img_spatial)
                            o = align_to_input(o, img_spatial);
                        end
                        V(:,:,:,1) = V(:,:,:,1) + o;   % reshape(o) broadcasts over [H,W,C]
                    end
                end
    
                % return output set
                image = ImageStar(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            else % reachability with Star sets
                % In general, this layer is used as the bias addition following a FullyConnectedLayer
                % Scale (weights)
                image = in_image; % create copy to perform operations (if need to)
                if obj.DoScale
                    % R2025b's ONNX importer emits ScalingLayer (mapped to
                    % ElementwiseAffineLayer) with Scale as a per-channel
                    % vector, not a square matrix. Wrap accordingly so Star
                    % dim is preserved (mirrors the ImageStar branch above).
                    % double() casts guard against single-precision ONNX scales.
                    s = obj.Scale;
                    if isscalar(s)
                        image = image.affineMap(double(s) * eye(image.dim), []);
                    elseif isvector(s)
                        image = image.affineMap(diag(double(s(:))), []);
                    else
                        image = image.affineMap(s, []); % W*x (full matrix)
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
                        % flatten produced k*dim copies. Taking the first dim
                        % slice is only correct if every tile is IDENTICAL --
                        % otherwise we'd silently use a wrong bias [13]. Verify
                        % the tiles match; refuse if they don't.
                        flat = double(o(:));
                        tiles = reshape(flat, image.dim, []);
                        if max(abs(tiles - tiles(:,1)), [], 'all') > 1e-9
                            error('ElementwiseAffineLayer:biasTileMismatch', ...
                                ['Offset has %d elements (%dx the Star dim %d) but the ' ...
                                 'tiles are not identical -- cannot safely reduce to a ' ...
                                 '[%d,1] bias without silently using a wrong value.'], ...
                                numel(o), numel(o)/image.dim, image.dim, image.dim);
                        end
                        o = tiles(:,1);
                    elseif numel(o) < image.dim && mod(image.dim, numel(o)) == 0
                        % A bias shorter than the Star dim that divides it is
                        % AMBIGUOUS [13/165]: repmat assumes a channel-fastest
                        % flat layout ([o;o;...]), but a position-fastest flatten
                        % needs repelem -- the layer has no tensor shape to decide.
                        % evaluate() itself errors on this flat case (a [C] bias
                        % cannot broadcast against a [dim] vector), so silently
                        % guessing a tiling makes reach disagree with evaluate.
                        % Refuse rather than guess (a per-channel bias on a spatial
                        % set goes through the shape-aware ImageStar branch above).
                        error('ElementwiseAffineLayer:biasLayoutAmbiguous', ...
                            ['Offset length %d divides Star dim %d but the tiling ' ...
                             'layout (channel- vs position-fastest) is ambiguous and ' ...
                             'evaluate() cannot broadcast it on a flat vector. Refusing ' ...
                             'to guess; supply a [dim,1] bias or use a shaped ImageStar.'], ...
                            numel(o), image.dim);
                    else
                        % numel(o) is neither the Star dim, a clean multiple, nor
                        % a clean divisor: padding with zeros or truncating would
                        % silently fabricate a wrong bias (unsound) [13]. Refuse.
                        error('ElementwiseAffineLayer:biasDimMismatch', ...
                            ['Offset length %d is incompatible with Star dim %d ' ...
                             '(not equal, multiple, or divisor). Refusing to pad/truncate ' ...
                             'into a wrong bias -- check the ONNX bias shape / flatten order.'], ...
                            numel(o), image.dim);
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
        function L = parse(src)
            % @src: either nnet.onnx.layer.ElementwiseAffineLayer (R2024b and
            %       earlier ONNX imports) or nnet.cnn.layer.ScalingLayer
            %       (R2025a+ ONNX imports -- semantically identical: Y = Scale.*X + Offset)
            % @L : a ElementwiseAffineLayer obj for reachability analysis purpose

            isElementwise = isa(src, 'nnet.onnx.layer.ElementwiseAffineLayer');
            isScaling     = isa(src, 'nnet.cnn.layer.ScalingLayer');
            if ~isElementwise && ~isScaling
                error(['Input must be nnet.onnx.layer.ElementwiseAffineLayer or ' ...
                       'nnet.cnn.layer.ScalingLayer (got %s).'], class(src));
            end

            % Scale / Offset are public on both classes and have the same semantics.
            scale  = src.Scale;
            offset = src.Offset;
            % ScalingLayer has no DoScale/DoOffset flags -- synthesize defaults
            % (the reach code treats both-true as "apply scale and offset", which
            % matches the ScalingLayer forward pass).
            if isprop(src, 'DoScale'),  doScale  = src.DoScale;  else, doScale  = true; end
            if isprop(src, 'DoOffset'), doOffset = src.DoOffset; else, doOffset = true; end

            L = ElementwiseAffineLayer(src.Name, scale, offset, doScale, doOffset);

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

function tf = elementwiseaffine_broadcastable(sz_p, sz_x)
% True when a parameter of size sz_p applies element-wise to an input of size
% sz_x WITHOUT changing x's size -- i.e. every parameter dim either equals the
% input dim or is a singleton that broadcasts over it. Same dimensionality is
% required so it applies as-is, without length-based re-alignment that could
% land its non-singleton dim on the wrong axis when two input dims share a size
% (e.g. a [1,1,C] channel bias on an [H,W,C] map with H==C). A flat/lower-rank
% parameter (e.g. [1,1,1,N] vs an [N,1] vector) has differing ndims -> false ->
% the caller falls back to align_to_input.
%
% CRITICAL: do NOT allow sz_x==1 with sz_p>1 -- that is broadcast-compatible in
% the MATLAB sense but EXPANDS x (e.g. a [1,5] row param on a [5,1] column input
% yields a [5,5] OUTER PRODUCT, not a [5,1] element-wise scale), which corrupts
% the downstream shape (the ACAS-Xu per-feature input-normalisation regression).
    if numel(sz_p) ~= numel(sz_x), tf = false; return; end
    tf = all(sz_p == sz_x | sz_p == 1);
end

