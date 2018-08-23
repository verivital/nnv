classdef FFNN
    % FeedForward Neural Network Class
    % Dung Tran: 8/22/2018
    
    properties
        Layers = []; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        nL = 0; % number of Layers
        nN = 0; % number of Neurons
        nI = 0; % number of Inputs
        nO = 0; % number of Outputs
    end
    
    methods
        
        % constructor
        function obj = FFNN(Layers)
            nL = size(Layers, 2); % number of Layer
            for i=1:nL
                L = Layers(i);
                if ~isa(L, 'Layer')
                    error('Element %d of Layers array is not a Layer object', i);
                end
            end
            
            % check consistency between layers
            for i=1:nL-1
                if size(Layers(i).W, 1) ~= size(Layers(i + 1).W, 2)
                    error('Inconsistent dimensions between Layer %d and Layer %d', i, i + 1);
                end
            end
            
            obj.Layers = Layers;
            obj.nL = nL;    % number of layers
            obj.nI = size(Layers(1).W, 2); % number of inputs
            obj.nO = size(Layers(nL).W, 1); % number of outputs
            
            nN = 0;
            for i=1:nL
                nN = nN + Layers(i).N;
            end
            obj.nN = nN; % number of neurons
            
        end
        
        % Evaluation of a FFNN
        function y = evaluate(obj, x)
            % Evaluation of this FFNN
            % @x: input vector x
            % @y: output vector y
            
            y = x; 
            for i=1:obj.nL
                y = obj.Layers(i).evaluate(y);
            end
        
        end
        
        % Reach Set computation
        function [R, rn, t] = reach(obj, I, option)
            % Reach Set Computation of this FFNN
            % @I: input set which is a polyhedron
            % @option: = 'exact' -> compute the exact reach set
            %          = 'approx' -> compute the over-approximate reach set
            % @R: output set which is an array of polyhedron if option = 'exact'
            %     or is a polyhedron if option = 'approx'
            % @rn: number of ReLU_i (stepReach) operation reduced
            % @t : computation time which is an array containing computation time
            %     on each layer.
        
            rn = [];
            t = [];
            
            R1 = I;
            for i=1:obj.nL
                [R1, rn1, t1] = obj.Layers(i).reach(R1, option);
                rn = [rn, rn1];
                t = [t, t1];
                fprintf('\nFinish reach set computation for layer %d in %.5f seconds', i, t1);
            end
            R = R1;
            
        end
    end
    
end

