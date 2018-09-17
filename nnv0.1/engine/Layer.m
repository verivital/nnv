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
            
            t1 = tic;
            
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
            
            t = toc(t1);              
        end
        
        % reach set computation in horizonal manner, contruct a box for the
        % output set of the layer
        
        function [R, t] = reach_approx_box(obj, I, nB, parallel)
            % @I : input set, a polyhedron
            % @nB: number of boxes used to over-approximate reach set
            % @parallel: @parallel = 'parallel' -> compute in parallel with multiple cores or
            % clusters
            %            @parallel = 'single' -> single core
            % @R: box
            % @t: computation time
            
            t1 = tic;
            n = length(I);
            R = [];
            if strcmp(obj.f, 'Linear')
               for i=1:n
                   % if activation function is linear we don't use box but
                   % return the exact output
                   R = [R I(i).affineMap(obj.W) + obj.b]; 
               end
                
                
            elseif strcmp(obj.f, 'ReLU')
                
                [R, rn, t2] = obj.reach_exact(I, parallel);
                R = Reduction.merge_box(R, nB, 'single');
                
            else                
                error('\nUnsupported activation function');
            end 
                
            t = toc(t1);
            
        end
        
        % over-approximate reach set using polyhedron without parallel
        % computing
        function [R, runtime] = reach_approx_polyhedron(obj, I, nP, nC, parallel)
            % @I : input set, a polyhedron
            % @R: over-approximate reach set (a polyhedron)
            % @t: computation time
            % @nP: number of polyhedra in the output R ( can set = number of cores)
            % @nC: number of constraints in each polyhedron of the output R
            %      @nC: is used to control the number of constraints of the
            %      polyhedra of the output, make it always ~nC.
            % @parallel: = 'parallel' -> using parallel computing
            %            = 'single' -> using single core for computing
            
            startTime = tic;
            
            for i=1:length(I)
                if I(i).isEmptySet
                    error('Input %d is an empty set', i);
                end
                if size(obj.W, 2) ~= size(I(i).A, 2)
                    error('Inconsistent dimensions between input set %d and weights matrix', i);
                end
                if ~I(i).isBounded
                    error('Input set is not bounded')
                end
            end
            
            if strcmp(obj.f, 'Linear')
                R = [];
                for i=1:length(I)
                    R = [R I(i).affineMap(obj.W) + obj.b];
                end                
            elseif strcmp(obj.f, 'ReLU')                
                
                % compute the exact reach set
                fprintf('\nComputing the exact reach set of the current layer...')
                [R, rn, t] = obj.reach_exact(I, parallel);
                fprintf('\nReach set computation is done in %.5f seconds', t);
                
                
                % clustering output reach set into nP groups based on their overlapness
                
                fprintf('\nClustering reach set (%d polyhedra) into %d groups based on their overlapness...', length(R), nP);
                t1 = tic;
                clustered_set = Reduction.clusterByOverlapness(R, nP); 
                t2 = toc(t1);
                fprintf('\nFinish clustering in %.5f seconds', t2);
                display(clustered_set);
                
                t1 = tic; 
                R = [];
                for j=1:nP
                    fprintf('\nMerging all polyhedra in the %dth clustered group into one polyhdron...', j);
                    % merge all polyhedra in one group into one polyhedron
                    R = [R Reduction.recursiveMerge(clustered_set{j, 1}, 1, parallel)]; 
                end
                t2 = toc(tic);
                fprintf('\nFinish merging in %.5f seconds', t2);
                
                % reduce the number of constraints of the output polyhedra
                t1 = tic;
                fprintf('\nReducing the constraints of the %d merged polyhedra to ~%d constraints', nP, nC)
                R = Reduction.reduceConstraints(R, nC, parallel); 
                t2 = toc(t1);
                fprintf('\nFinish constraints reduction in %.5f seconds', t2);
            else
                error('Unsupported activiation function');
            end
            
            runtime = toc(startTime);                                          
            
        end
        
        
        % coarse over-approximation of reachable set using box
        
        function [R, runtime] = reach_approx_box_coarse(obj, I, parallel)
            % @I: input array
            % @parallel: = 'parallel' using parallel computing
            %            = 'single' using single core for computing
            
            t1 = tic;
            n = length(I);
            R = [];
            
            if strcmp(parallel, 'single')
                if strcmp(obj.f, 'Linear')
                    for i=1:n
                        R = [R I(i).affineMap(obj.W) + obj.b];
                    end
                elseif strcmp(obj.f, 'ReLU')
                    
                    for i=1:n
                        I1 = I(i).affineMap(obj.W) + obj.b;
                        I1.outerApprox;
                        lb = I1.Internal.lb;
                        ub = I1.Internal.ub;

                        for j=1:length(lb)
                            if lb(j) < 0
                                lb(j) = 0;
                            end
                            if ub(j) < 0
                                ub(j) = 0;
                            end
                        end

                        R = [R Polyhedron('lb', lb, 'ub', ub)];

                    end
                else
                    error('Unsupported activation function');
                end
                
            
            elseif strcmp(parallel, 'parallel')
                
                
                if strcmp(obj.f, 'Linear')
                    parfor i=1:n
                        R = [R I(i).affineMap(obj.W) + obj.b];
                    end
                elseif strcmp(obj.f, 'ReLU')
                    
                    parfor i=1:n
                        I1 = I(i).affineMap(obj.W) + obj.b;
                        I1.outerApprox;
                        lb = I1.Internal.lb;
                        ub = I1.Internal.ub;

                        for j=1:length(lb)
                            if lb(j) < 0
                                lb(j) = 0;
                            end
                            if ub(j) < 0
                                ub(j) = 0;
                            end
                        end

                        R = [R Polyhedron('lb', lb, 'ub', ub)];

                    end
                    
                else
                    error('Unsupported activation function');
                end
                
            else 
                error('Unknown parallel computing option');
            end
                
            runtime = toc(t1);
            
            
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

