classdef Layer
    % A base class for a layer in NN
    % Dung Tran: 8/20/2018
    
    properties
        W; % weights_mat 
        b; % bias vector 
        f; % activation function;
        N; % number of neurons
    end
    
    methods
        % Constructor
        function obj = Layer(W, b, f)
            if size(W, 1) == size(b, 1)
                obj.W = W;
                obj.b = b;
                obj.N = size(W, 1);
            else
                error('Inconsistent dimensions between Weights matrix and bias vector');
            end
                
            if ~ischar(f) && ~strcmp(f, 'ReLU') && ~strcmp(f, 'Linear') && ~strcmp(f, 'Tanh')
                error('Invalid or Unsupported activation function');
            else
                obj.f = f;
            end
               
        end
        
        % Evaluation method
        function y = evaluate(obj, x)  % evaluation of this layer with a specific vector
            if size(x, 1) ~= size(obj.W, 2) || size(x, 2) ~= 1
                error('Invalid or inconsistent input vector')
            end
            y1 = obj.W * x + obj.b;
            y = zeros(obj.N, 1);
            for i = 1:obj.N
                if strcmp(obj.f, 'ReLU')
                    if y1(i, 1) >= 0
                        y(i, 1) = y1(i, 1);
                    else
                        y(i, 1) = 0;
                    end
                elseif strcmp(obj.f, 'Linear')
                    y(i, 1) = y1(i, 1);
                elseif strcmp(obj.f, 'Tanh')
                    y(i, 1) = tanh(y1(i, 1));
                end 
            end
        end
        
        % Reachable Set Computation
        function [R, rn, t] = reach(obj, inputSetArray, option)
            % compute reachable set of a layer
            % @inputSetArray: an array of input sets which are bounded polyhedra
            % @option: = 'exact' or 'approx'
            % if option = exact: exact reachable set is computed
            % if option = 'approx': over-approximate reachable set is
            % computed.
            % @R: the reachable set
            % @rn: the number of ReLU operation reduced
            % @t : Elapsed time for reach set computation
            
            tic;
            
            rn = 0; 
            R = [];
            p = length(inputSetArray); % number of polyhedron in the input set
            for i=1:p
                I = inputSetArray(i); 
                if size(obj.W, 2) ~= size(I.A, 2)
                    error('Inconsistent dimensions between input set and weights matrix')
                end

                I1 = I.affineMap(obj.W); % affine map Wx, x in I
                I1 = I1 + obj.b; % affine map Wx + b, x in I

                if strcmp(obj.f, 'Linear')
                    R1 = I1;
                    rn1 = 0;
                elseif strcmp(obj.f, 'ReLU')
                    [R1, rn1] = ReLU.reach(I1, option);
                else
                    error('Unsupported activation function, currently support ReLU and Linear')
                end
                R = [R, R1];
                rn = rn + rn1;
            end            
            t = toc;              
        end
        
        % evaluate the value of the layer output with all vertices of input set 
        function Y = sample(obj, V)
  
            n = size(V, 2);
            Y = [];
            for j=1:n
                y = obj.evaluate(V(:, j));
                Y = [Y y];
            end
        end
        
               
    end
    
end

