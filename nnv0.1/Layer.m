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
                
            if ~ischar(f) || ~strcmp(f, 'ReLU') || ~strcmp(f, 'Linear') || ~strcmp(f, 'Tanh')
                error('Invalid activation function');
            else
            obj.f = f;
            end
               
        end
        
        % Evaluation method
        function y = evaluation(obj, x)  % evaluation of this layer with a specific vector
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
        
        
    end
    
end

