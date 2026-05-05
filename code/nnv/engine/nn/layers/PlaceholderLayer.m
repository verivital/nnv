classdef PlaceholderLayer < handle
    % Placeholder Layer object
    % This layer serves as a placeholder for all the layers that to not
    % have any effect on the verification algorithm such as Dropout or
    % RegressionLayer from MATLAB
    %
    % Author: Diego Manzanas Lopez
    % Date: 03/15/2023
    
    properties
        Name = 'NoOpLayer';
        Type = ''; % e.g. dropout
        Perm = [];  % optional MATLAB-style permutation order (1-indexed)
                    % set when this placeholder represents a Transpose op
    end

    methods % constructor

        % create layer
        function obj = PlaceholderLayer(name, Ltype)
            % @name: name of the layer
            % @type: original layer type from MATLAB
            obj.Name = name;
            obj.Type = Ltype;
        end

    end

    methods % main methods

        % evaluate
        function out_im = evaluate(obj, inputs)
            % return output = input, optionally permuted or with an
            % element-wise function applied (Sign/Abs/etc.).
            if ~isempty(obj.Perm)
                % MATLAB's permute requires perm to cover at least ndims(inputs).
                % If the perm is shorter (e.g. ONNX rank-2 perm but MATLAB
                % already promoted the tensor to rank-3), pad with identity.
                p = obj.Perm(:).';
                k = max(numel(p), ndims(inputs));
                if numel(p) < k
                    p = [p, (numel(p)+1):k];
                end
                out_im = permute(inputs, p);
                return;
            end
            % Element-wise op stored as the Type tag (set by the loader for
            % ops like Sign / Abs that NNV doesn't have a dedicated layer for).
            switch obj.Type
                case 'Sign'
                    out_im = sign(inputs);
                case 'Abs'
                    out_im = abs(inputs);
                otherwise
                    out_im = inputs;
            end
        end

        function out_sq = evaluateSequence(~, inputs)
            out_sq = inputs;
        end
        
        % reachability analysis with multiple inputs
        function IS = reach(varargin)
            % return input set as output set
            % Regardless of options, input to layer corresponds to varargin{2}
            in_images = varargin{2}; 
            IS = in_images;
        end

        % reachability analysis with multiple inputs
        function IS = reachSequence(varargin)
            % return input set as output set
            % Regardless of options, input to layer corresponds to varargin{2}
            in_seqs = varargin{2}; 
            IS = in_seqs;
        end

    end

    methods % helper method

        % change params to gpuArrays
        function obj = toGPU(obj)
            % nothing to change in here (no params)
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, ~)
            % nothing to change in here (no params)
        end
        
    end
    
    methods(Static)
        
        % parsing method
        function L = parse(layer)
            % create NNV layer 
            L = PlaceholderLayer(layer.Name, class(layer));
        end

    end

end



