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
        
        % Exact reachable Set computation
        function [R, rn, t] = reach_exact(obj, inputSetArray, parallel)
            % compute reachable set of a layer
            % @inputSetArray: an array of input sets which are bounded polyhedra
            % @option: = 'exact' or 'approx'
            % computed.
            % @R: the reachable set
            % @rn: the number of ReLU operation reduced
            % @t : Elapsed time for reach set computation
            % @parallel: = 'parallel' use parallel computing
            %            = 'single' use single core
            
            tic;
            
            rn = 0; 
            R = [];
            p = length(inputSetArray); % number of polyhedron in the input set
            
            if strcmp(parallel, 'parallel') % use parallel computing                
            
                parfor i=1:p
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
                        [R1, rn1] = ReLU.reach(I1);

                    else
                        error('Unsupported activation function, currently support ReLU and Linear')
                    end
                    R = [R, R1];
                    rn = rn + rn1;
                end
                
            elseif strcmp(parallel, 'single') % use single core for computing
                
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
                        [R1, rn1] = ReLU.reach(I1);

                    else
                        error('Unsupported activation function, currently support ReLU and Linear')
                    end
                    R = [R, R1];
                    rn = rn + rn1;
                end
                
            else                
                error('Unknown computation option: parallel or single');
            end
            
            t = toc;              
        end
        
        % reach set computation in horizonal manner, contruct a box for the
        % output set of the layer
        
        function [R, t] = reach_approx_box(obj, I, parallel)
            % @I : input set, a polyhedron
            % @parallel: @parallel = 'parallel' -> compute in parallel with multiple cores or
            % clusters
            %            @parallel = 'single' -> single core
            % @R: box
            % @t: computation time
            
            tic;
            
            if size(obj.W, 2) ~= size(I.A, 2)
                    error('Inconsistent dimensions between input set and weights matrix')
            end
                                    
            R = I.affineMap(obj.W) + obj.b;
            
            if strcmp(obj.f, 'Linear')
               % do nothing 
            elseif strcmp(obj.f, 'ReLU')
                
                R.outerApprox;
                lb = zeros(obj.N, 1);
                ub = zeros(obj.N, 1);
                
                if strcmp(parallel, 'parallel') % computing in parallel
                                    
                    lb1 = R.Internal.lb;
                    ub1 = R.Internal.ub;
                    
                    N = obj.N; % this needs to be done for parallel computing                   
                    parfor i=1:N
                        if lb1(i) <= 0
                            lb(i) = 0;
                        else
                            lb(i) = lb1(i);
                        end
                        if ub1(i) <= 0
                            ub(i) = 0;
                        else
                            ub(i) = ub1(i);
                        end
                    end
                    
                elseif strcmp(parallel, 'single')
                    
                    lb = R.Internal.lb;
                    ub = R.Internal.ub;
                    
                    for i=1:obj.N %single core computing
                        if lb(i) <= 0
                            lb(i) = 0;
                        end
                        if ub(i) <= 0
                            ub(i) = 0;
                        end
                    end
                    
                else 
                    error('\nUnknown option');                  
                end
                
                
                R = Polyhedron('lb', lb, 'ub', ub);
                
            else
                
                error('\nUnsupported activation function');
                
            end 
                
            t = toc;
            
        end
        
        % over-approximate reach set using polyhedron without parallel
        % computing
        function [R, runtime] = reach_approx_polyhedron(obj, I, parallel)
            % @I : input set, a polyhedron
            % @R: over-approximate reach set (a polyhedron)
            % @t: computation time
            % @parallel: = 'parallel' -> using parallel computing
            %            = 'single' -> using single core for computing
            
            startTime = tic;
            
            if I.isEmptySet
                error('Input set is an empty set');
            end
            
            if size(obj.W, 2) ~= size(I.A, 2)
                    error('Inconsistent dimensions between input set and weights matrix')
            end
           
            
            if ~I.isBounded
                error('Input set is not bounded')
            end
            
            if strcmp(obj.f, 'Linear')
                R = I.affineMap(obj.W) + obj.b;
            elseif strcmp(obj.f, 'ReLU')                
                R = ReLU.reach(I.affineMap(obj.W) + obj.b);
                R = Reduction.recursiveMerge(R, parallel);                                
            else
                error('Unsupported activiation function');
            end
            
            runtime = toc(startTime);                                          
            
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
    
    methods(Static) % Static methods for cluster computing
        
        
    end
    
end

