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
        function out_im = evaluate(~, inputs)
            % return output = input
            out_im = inputs;
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
    
    methods(Static)
        
        % parsing method
        function L = parse(layer)
            % create NNV layer 
            L = PlaceholderLayer(layer.Name, class(layer));
        end

    end

end



