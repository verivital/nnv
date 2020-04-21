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
                obj.W = double(W);
                obj.b = double(b);
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
                display(x);
                display(obj.W);
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
        
        % Exact reachable Set computation for one layer using Polyhedron or
        % Star
        function [R, rn, t] = reach_exact(obj, inputSetArray, parallel)
            % @inputSetArray: an array of input sets which are bounded polyhedra
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
                
                if p > 1
                    parfor i=1:p
                        I = inputSetArray(i); 
                        if isa(I, 'Polyhedron')                            
                            I1 = I.affineMap(obj.W, 'vrep') + obj.b;
                        elseif isa(I, 'Star')
                            I1 = I.affineMap(obj.W, obj.b);
                        else
                            error('Unknown set type');
                        end
                      
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
                    
                else  % if single input, use ReLU parallel
                    
                    I = inputSetArray(1);
                    if isa(I, 'Polyhedron')                            
                        I1 = I.affineMap(obj.W, 'vrep') + obj.b;
                    elseif isa(I, 'Star')
                        I1 = I.affineMap(obj.W, obj.b);
                    else
                        error('Unknown set type');
                    end
           
                    if strcmp(obj.f, 'Linear')
                        R1 = I1;
                        rn1 = 0;
                    elseif strcmp(obj.f, 'ReLU')
                        [R1, rn1] = ReLU.reach_parallel(I1);
                    else
                        error('Unsupported activation function, currently support ReLU and Linear')
                    end
                    R = [R, R1];
                    rn = rn + rn1;
                    
                end
                
                
            elseif strcmp(parallel, 'single') % use single core for computing
                
                for i=1:p
                    I = inputSetArray(i);
   
                    if isa(I, 'Polyhedron')                            
                        I1 = I.affineMap(obj.W, 'vrep') + obj.b;
                    elseif isa(I, 'Star')
                        I1 = I.affineMap(obj.W, obj.b);
                    else
                        error('Unknown set type');
                    end
                    
                    if strcmp(obj.f, 'Linear')
                        R1 = I1;
                        rn1 = 0;
                    elseif strcmp(obj.f, 'ReLU')
                        [R1, rn1] = ReLU.reach(I1);
                        %[R1, rn1] = ReLU.reach_new(I1);

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
                
        
        % coarse over-approximation of reachable set using box
        
        function [R, runtime] = reach_approx_box_coarse(obj, I, nP, parallel)
            % @I: input array
            % @nP: maximum number of polyhedra to over-approximate reach
            % set for each layer
            % @parallel: = 'parallel' using parallel computing
            %            = 'single' using single core for computing
            
            t1 = tic;
            n = length(I);
            R = [];
            
            if strcmp(parallel, 'single')
                if strcmp(obj.f, 'Linear')
                    
                    if ~isa(I(1), 'Box')
                        for i=1:length(I)
                            R = [R I(i).affineMap(obj.W) + obj.b];
                        end
                        
                    else
                        for i=1:length(I)
                            R = [R I(i).affineMap(obj.W, obj.b)];
                        end
                        
                    end
                                        
                elseif strcmp(obj.f, 'ReLU')
                    
                    for i=1:n
                        
                        if ~isa(I(i), 'Box')
                            I1 = Reduction.affineMap(I(i), obj.W) + obj.b;
                            I1.outerApprox;
                            lb = I1.Internal.lb;
                            ub = I1.Internal.ub;
                            
                        else 
                            
                            I1 = I(i).affineMap(obj.W, obj.b);
                            lb = I1.lb;
                            ub = I1.ub;
                       
                        end

                        for j=1:length(lb)
                            if lb(j) < 0
                                lb(j) = 0;
                            end
                            if ub(j) < 0
                                ub(j) = 0;
                            end
                        end

                        R = [R Box(lb, ub)];     
       
                    end
                    
                    if ~isempty(nP) && length(R) > nP
                        R = Reduction.merge_box2(R, nP, parallel);
                    end
                else
                    error('Unsupported activation function');
                end
                
            
            elseif strcmp(parallel, 'parallel')
                
                
                if strcmp(obj.f, 'Linear')
                    
                    if ~isa(I(1), 'Box')
                        for i=1:length(I)
                            R = [R I(i).affineMap(obj.W) + obj.b];
                        end
                        
                    else
                        for i=1:length(I)
                            R = [R I(i).affineMap(obj.W, obj.b)];
                        end
                        
                    end
                    
                elseif strcmp(obj.f, 'ReLU')
                    
                    parfor i=1:n
                        
                        if ~isa(I(i), 'Box')
                            
                            I1 = Reduction.affineMap(I(i), obj.W) + obj.b;
                            I1.outerApprox;
                            lb = I1.Internal.lb;
                            ub = I1.Internal.ub;
                            
                        else 
                            
                            I1 = I(i).affineMap(obj.W, obj.b);
                            lb = I1.lb;
                            ub = I1.ub;
                        end

                        for j=1:length(lb)
                            if lb(j) < 0
                                lb(j) = 0;
                            end
                            if ub(j) < 0
                                ub(j) = 0;
                            end
                        end

                        R = [R Box(lb, ub)];     

       
                    end
                    
                    if ~isempty(nP) && length(R) > nP
                        R = Reduction.merge_box2(R, nP, parallel);
                    end
                    
                else
                    error('Unsupported activation function');
                end
                
            else 
                error('Unknown parallel computing option');
            end
                
            runtime = toc(t1);
            
            
        end
        
        
        % reach set computation by mixing between boxes and polyhedra
        
        function [P, runtime] = reach_mix(obj, I, nP, parallel)
            % @I: input polyhedra
            % @nP: maximum allowable total number of polyhedra at output
          
            
            % @parallel: = 'parallel' using parallel computing
            %            = 'single' using single core computing
            
           t1 = tic; 
           n = length(I);
           
           if length(I) < nP                        
               P = obj.reach_exact(I, parallel);
           else
               if ~isa(I(1), 'Box')
                   P = Reduction.merge_box(I, nP, parallel);
               else
                   P = Reduction.merge_box2(I, nP, parallel);
               end
               
               P = obj.reach_approx_box_coarse(P, nP, parallel);
           end
           
           runtime = toc(t1);
            
            
        end
        
        
        % reachable set computation with controlling number of stars
        function [P, rn, runtime] = reach_approx_star(obj, I, parallel)
            
            % @I: input star
            % @parallel: = 'parallel' using parallel computing
            %            = 'single' using single core computing
            
            % author: Dung Tran
            % date: 2/25/2019
            
            t1 = tic; 
            n = length(I);
           
            if ~isa(I, 'Star')
                error('Input is not a Star');
            else
                
                I1 = I.affineMap(obj.W, obj.b);
                if strcmp(obj.f, 'Linear')
                    P = I1;
                    rn = 0;
                elseif strcmp(obj.f, 'ReLU')
                    [P, rn] = ReLU.reach_approx(I1, parallel);
                else
                    error('Unsupported activation function, currently support ReLU and Linear')
                end
                                
                
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

