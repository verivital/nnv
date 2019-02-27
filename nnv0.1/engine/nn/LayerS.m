classdef LayerS
    % LAYERS is a class that contains only reachability analysis method using
    % stars, this gonna replace Layer class in the future, i.e., delete
    % polyhedron-based reachability analysis method.
    % author: Dung Tran
    % date: 27/2/2019
    
    properties
        W; % weights_mat 
        b; % bias vector 
        f; % activation function;
        N; % number of neurons
    end
    
    methods % constructor - evaluation - sampling
        % Constructor
        function obj = LayerS(W, b, f)
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
                elseif strcmp(obj.f, 'SatLin')
                    y(i, 1) = satlin(y1(i, 1));             
                end 
            end
        end
        
        % evaluate the value of the layer output with a set of vertices
        function Y = sample(obj, V)
  
            n = size(V, 2);
            Y = [];
            for j=1:n
                y = obj.evaluate(V(:, j));
                Y = [Y y];
            end
        end
        
    end
    
    
    methods % reachability analysis method
        
        function S = reach(varargin)
            % @I: an array of inputs
            % @method: 'star' or 'abs-dom', i.e., abstract domain (support
            % later) or 'face-latice', i.e., face latice (support later)
            % @option:  'parallel' use parallel computing
            %           '[]' or not declared -> don't use parallel
            %           computing
            
            % author: Dung Tran
            % date: 27/2/2019
            
             
            % parse inputs 
            switch nargin
                case 4
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                case 3
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 3 or 4)');
            end
            
            S = [];
            
            % do reachability analysis
            if strcmp(method, 'star') % reachability analysis using star set
                
                n = length(I);
                for i=1:n
                    
                    % affine mapping y = Wx + b;
                    if isa(I(i), 'Star')
                        I1 = I(i).affineMap(obj.W, obj.b);
                    else
                        error('%d^th input is not a star', i);
                    end
                    % apply activation function: y' = ReLU(y) or y' = 
                    
                    if strcmp(obj.f, 'purelin')
                        S = [S I1];
                    elseif strcmp(obj.f, 'poslin')
                        S = [S PosLin.reach(I1, option)];
                     elseif strcmp(obj.f, 'satlin')
                        S = [S SatLin.reach(I1, option)];
                    else
                        error('Unsupported activation function, currently support pureline, posline, satlin')
                    end
                    
                    
                end
                
                
            elseif strcmp(method, 'abs-dom') % reachability analysis using abstract-domain
                
                error('Abstract-domain method is not supported yet');
                
            elseif strcmp(method, 'face-latice') % reachability analysis using face-latice
                
                error('Face-latice method is not supported yet');
                
            else
                error('Unknown method');
            end
            
        end
    
    end
    
    
    
end

