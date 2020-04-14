classdef PosLin
    % POSLIN Class contains method for reachability analysis for Layer with ReLU activation function 
    % IT GONNA BE USED TO REPLACE ReLU class in the future
    % Its name is consistent with the name of ReLU acivation function in
    % matlab called 'poslin'
    
    % author: Dung Tran
    % date: 27/2/2019
    
    properties
        
    end
    
    methods(Static) % evaluate method and reachability analysis with stars 
        
        % evaluation
        function y = evaluate(x)
            y = poslin(x);
        end
        
        % stepReach method, compute reachable set for a single step
        function S = stepReach(I, index)
            % @I: single star set input
            % @index: index of the neuron performing stepPosLin
            % @xmin: minimum of x[index]
            % @xmax: maximum of x[index]
            % @S: star output set
            
            % author: Dung Tran
            % date: 27/2/2019
            
            if ~isa(I, 'Star')
                error('Input is not a star set');
            end
            
            [xmin, xmax] = I.getRange(index);
                       
            if xmin >= 0
                S = I; 
            elseif xmax <= 0 
                
                V1 = I.V;
                V1(index, :) = 0;
                if ~isempty(I.Z)
                    c = I.Z.c;
                    c(index) = 0;
                    V = I.Z.V;
                    V(index, :) = 0;
                    new_Z = Zono(c, V); % update outer-zono
                else
                    new_Z = [];
                end
                S = Star(V1, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);
                
            elseif xmin < 0 && xmax > 0
                % S1 = I && x[index] < 0 
                c = I.V(index, 1);
                V = I.V(index, 2:I.nVar + 1); 
                new_C = vertcat(I.C, V);
                new_d = vertcat(I.d, -c);                
                new_V = I.V;
                new_V(index, :) = zeros(1, I.nVar + 1);
                
                % update outer-zono
                if ~isempty(I.Z)
                    c1 = I.Z.c;
                    c1(index) = 0;
                    V1 = I.Z.V;
                    V1(index, :) = 0;
                    new_Z = Zono(c1, V1);
                else
                    new_Z = [];
                end
                S1 = Star(new_V, new_C, new_d, I.predicate_lb, I.predicate_ub, new_Z);
                
                % S2 = I && x[index] >= 0
                new_C = vertcat(I.C, -V);
                new_d = vertcat(I.d, c);
                
                S2 = Star(I.V, new_C, new_d, I.predicate_lb, I.predicate_ub, I.Z);
                
%                 a = S1.isEmptySet;
%                 b = S2.isEmptySet;
%                                              
%                 if a && ~b
%                     S = S2;
%                 end
%                 if a && b
%                     S = [];
%                 end
%                 if ~a && b
%                     S = S1;
%                 end
%                 if ~a && ~b
%                     S = [S1 S2];
%                 end

               S = [S1 S2];
        
            end
                        
        end
        
        
        % stepReach with multiple inputs
        function S = stepReachMultipleInputs(I, index, option)
            % @I: an array of stars
            % @index: index where stepReach is performed
            % @option: = 'parallel' use parallel computing
            %          = not declare -> don't use parallel computing
            
            % author: Dung Tran
            % date: 27/2/2019
            % update: 4/2/2020
            
            p = length(I);
            S = [];
            
            if isempty(option)
                
                for i=1:p
                    S =[S, PosLin.stepReach(I(i), index)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    S =[S, PosLin.stepReach(I(i), index)];
                end
                
            else
                error('Unknown option');
            end
            
            
        end
        
        
        % exact reachability analysis using star
        function S = reach_star_exact(I, option)
            % @I: star input sets
            % @option: = 'parallel' using parallel computing
            %          = '[]'    do not use parallel computing
            
            % author: Dung Tran
            % date: 3/16/2019
            
            
             if ~isempty(I)
                            
                %[lb, ub] = I.getRanges;
                [lb, ub] = I.estimateRanges;
                
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    map = find(ub <= 0); % computation map
                    V = I.V;
                    V(map, :) = 0;
                    % update outer-zono
                    if ~isempty(I.Z)
                        c1 = I.Z.c;
                        c1(map, :) = 0;
                        V1 = I.Z.V;
                        V1(map, :) = 0;
                        new_Z = Zono(c1, V1);
                    else
                        new_Z = [];
                    end
                    
                    In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);                    
                    map = find(lb <0 & ub > 0);
                    m = length(map);                    
                    for i=1:m
                        fprintf('\nPerforming exact PosLin_%d operation using Star', map(i));
                        In = PosLin.stepReachMultipleInputs(In, map(i), option);
                    end               
                    S = In;
                end
                
            else
                S = [];
            end
            
        end
        
                % exact reachability analysis using star
        function S = reach_star_exact_multipleInputs(In, option)
            % @In: star input sets
            % @option: = 'parallel' using parallel computing
            %          = '[]'    do not use parallel computing
            
            % author: Dung Tran
            % date: 7/25/2019
            
            
             n = length(In);
             S = [];
             if strcmp(option, 'parallel')
                 parfor i=1:n
                     S = [S PosLin.reach_star_exact(In(i), [])];
                 end
             elseif isempty(option) || strcmp(option, 'single')
                 for i=1:n
                     S = [S PosLin.reach_star_exact(In(i), [])];
                 end
             else
                 error('Unknown computation option');
             end
        end
              
        
        
        % step reach approximation using star
        function S = stepReachStarApprox(I, index)
            % @I: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Dung Tran
            % date: 7/22/2019
                       
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end
            
            [lb, ub] = I.getRange(index);
           
            if lb >= 0
                S = I;
            elseif ub <= 0
                V = I.V;
                V(index, :) = 0;
                if ~isempty(I.Z)
                    c1= I.Z.c;
                    c1(index) = 0;
                    V1 = I.Z.V;
                    V1(index, :) = 0;
                    new_Z = Zono(c1, V1); % update outer-zono
                else
                    new_Z = [];
                end
                S = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);
            elseif lb < 0 && ub > 0
                fprintf('\nAdd a new predicate variables at index = %d', index);
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y[index] = ReLU(x[index]) >= 0
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                % constraint 2: y[index] >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                % constraint 3: y[index] <= ub(x[index] -lb)/(ub - lb)
                C3 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                d3 = -ub*lb/(ub-lb) +  ub*I.V(index, 1)/(ub-lb);

                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;              
                new_predicate_lb = [I.predicate_lb; 0];                
                new_predicate_ub = [I.predicate_ub; ub];
                
                % update outer-zono
                
                lamda = ub/(ub -lb);
                mu = -0.5*ub*lb/(ub - lb);
                if ~isempty(I.Z)
                    c = I.Z.c; 
                    c(index) = lamda * c(index) + mu;
                    V = I.Z.V;
                    V(index, :) = lamda * V(index, :);
                    I1 = zeros(I.dim,1);
                    I1(index) = mu;
                    V = [V I1];              
                    new_Z = Zono(c, V);
                else
                    new_Z = [];
                end
                
                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub, new_Z);
            end

        end
        
        
        % over-approximate reachability analysis using Star
        function S = reach_star_approx(I)
            % @I: star input set
            % @S: star output set

            % author: Dung Tran
            % date: 4/3/2019

            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isempty(I)
                S = [];
            else
                [lb, ub] = I.estimateRanges;
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    map = find(ub <= 0); % computation map
                    V = I.V;
                    V(map, :) = 0;
                    if ~isempty(I.Z)
                        c1 = I.Z.c;
                        c1(map, :) = 0;
                        V1 = I.Z.V;
                        V1(map, :) = 0;
                        new_Z = Zono(c1, V1);
                    else
                        new_Z = [];
                    end
                    In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);                    
                    map = find(lb <0 & ub > 0);
                    m = length(map); 
                    for i=1:m
                        fprintf('\nPerforming approximate PosLin_%d operation using Star', map(i));
                        In = PosLin.stepReachStarApprox(In, map(i));
                    end
                    S = In;
             
                end
            end

        end

    end
    
    
    methods(Static) % reachability analysis using Polyhedron
        
        % step reach set y = ReLU(x)
        function R = stepReach_Polyhedron(I, index, xmin, xmax)
            % @I : input set, a polyhedron
            % @i : index of current x[index] of current step
            % @xmin: min of x[index]
            % @xmax: max of x[index]
           
            I.normalize;
            dim = I.Dim;
            if xmin >= 0
                R = I; 
            elseif xmax < 0 
                Im = eye(dim);
                Im(index, index) = 0;
                R = I.affineMap(Im, 'vrep');
            elseif xmin < 0 && xmax >= 0
                
                Z1 = zeros(1, dim);
                Z1(1, index) = 1;
                Z2 = zeros(1, dim);
                Z2(1, index) = -1;

                A1 = vertcat(I.A, Z1);
                A2 = vertcat(I.A, Z2);
                b  = vertcat(I.b, [0]);
                R1 = Polyhedron('A', A1, 'b', b, 'Ae', I.Ae, 'be', I.be);
                R2 = Polyhedron('A', A2, 'b', b, 'Ae', I.Ae, 'be', I.be);
                
                Im = eye(dim);
                Im(index, index) = 0;
                R1 = R1.affineMap(Im, 'vrep');
                if R1.isEmptySet 
                    if R2.isEmptySet
                        R = [];
                    else
                        
                        R = R2;
                    end
                else
                    if R2.isEmptySet
                        R = R1;
                    else
                        R = [R1 R2];
                    end
                    
                end               
                
            end
            
        end
        
        
        % stepReach for multiple Input Sets 
        function R = stepReachMultipleInputs_Polyhedron(varargin)
            % @I: an array of input sets which are polyhedra
            % @index: index of current x[index] of current step
            % @xmin: min value of x[index]
            % @xmax: max value of x[index]
            % @option: = 'exact' -> compute an exact reach set
            %          = 'approx' -> compute an over-approximate reach set
            
            switch nargin
                
                case 5
                    I = varargin{1};
                    index = varargin{2};
                    xmin = varargin{3};
                    xmax = varargin{4};
                    option = varargin{5};
                
                otherwise
                    error('Invalid number of input arguments (should be 5)');
            end
            
            
            
            p = length(I);
            R = [];
            
            if isempty(option)
                
                for i=1:p
                    R =[R, PosLin.stepReach_Polyhedron(I(i), index, xmin, xmax)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    R =[R, PosLin.stepReach_Polyhedron(I(i), index, xmin, xmax)];
                end
                
            else
                error('Unknown option');
            end
            
            
                     
        end
        
        
        % exact reachability analysis using Polyhedron
        function R = reach_polyhedron_exact(I, option)
            
                        % @I: star input sets
            % @option: = 'parallel' using parallel computing
            %          = ''    do not use parallel computing
            
            % author: Dung Tran
            % date: 3/16/2019
            
             if ~isempty(I)
                if isa(I, 'Polyhedron')            
                    I.outerApprox; % find bounds of I state vector
                    lb = I.Internal.lb; % min-vec of x vector
                    ub = I.Internal.ub; % max-vec of x vector
                else
                    error('Input set is not a Polyhedron');
                end
                
                if isempty(lb) || isempty(ub)
                    R = [];
                else
                    map = find(lb < 0); % computation map
                    m = size(map, 1); % number of stepReach operations needs to be executed
                    In = I;
                    for i=1:m
                        fprintf('\nPerforming exact PosLin_%d operation using Polyhedron', map(i));
                        In = PosLin.stepReachMultipleInputs_Polyhedron(In, map(i), lb(map(i)), ub(map(i)), option);
                    end               
                    R = In;
                end
                
            else
                R = [];
            end
                 
        end
             
    end


    methods(Static) % over-approximate reachability analysis use zonotope
        
        % step over-approximate reachability analysis using zonotope
        function Z = stepReachZonoApprox(I, index, lb, ub)
            % @I: zonotope input set
            % @lb: lower bound of input at specific neuron i
            % @ub: lower bound of input at specfic neuron i
            % @index: index of the neuron we want to perform stepReach
            % @Z: zonotope output set
            
            % author: Dung Tran
            % date: 5/3/2019
            % update: 1/4/2020:
            %     update reason: change getBounds method for Zonotope
        
            % reference: Fast and Effective Robustness Ceritification,
            % Gagandeep Singh, NIPS 2018
           
            if ~isa(I, 'Zono')
                error('Input is not a Zonotope');
            end
            
            if lb >= 0
                Z = Zono(I.c, I.V);
                
            elseif ub <= 0
                c = I.c;
                c(index) = 0;
                V = I.V;
                V(index, :) = zeros(1, size(I.V, 2));
                Z = Zono(c, V);
                
            elseif lb < 0 && ub > 0
                
                lamda = ub/(ub -lb);
                mu = -0.5*ub*lb/(ub - lb);               
                
                c = I.c; 
                c(index) = lamda * c(index) + mu;
                V = I.V;
                V(index, :) = lamda * V(index, :);
                I1 = zeros(I.dim,1);
                I1(index) = mu;
                V = [V I1];
                
                Z = Zono(c, V);                
            end
            
        end
            
            
        % over-approximate reachability analysis use zonotope
        function Z = reach_zono_approx(I)
            % @I: zonotope input
            % @Z: zonotope output
            
            % author: Dung Tran
            % date: 5/3/2019
            
            % reference: Fast and Effective Robustness Ceritification,
            % Gagandeep Singh, NIPS 2018
            
            if ~isa(I, 'Zono')
                error('Input is not a Zonotope');
            end
                      
            In = I;
            [lb, ub] = I.getBounds;
            for i=1:I.dim
                fprintf('\nPerforming approximate PosLin_%d operation using Zonotope', i);
                In = PosLin.stepReachZonoApprox(In, i, lb(i), ub(i));
            end
            Z = In;
                       
        end
        
        
    end
    
    
    methods(Static) % reachability analysis using abstract-domain
        
        % step over-approximate reachability analysis using abstract-domain
        % we use star set to represent abstract-domain
        function S = stepReachAbstractDomain(varargin)
            % @I: star-input set
            % @index: index of neuron performing stepReach
            % @S: star output set represent abstract-domain of the output
            % set
            
            % author: Dung Tran
            % date: 5/3/2019
        
            % reference: An Abstract Domain for Certifying Neural Networks,
            % Gagandeep Singh, POPL 2019
            
            
            switch nargin
                
                case 4
                    
                    I = varargin{1};
                    index = varargin{2};
                    lb = varargin{3};
                    ub = varargin{4};
                
                case 2
                    I = varargin{1};
                    index = varargin{2};
                    %[lb, ub] = I.getRange(index); % our improved approach
                    [lb, ub] = I.estimateRange(index); % originial DeepPoly approach use estimated range
                otherwise
                    error('Invalid number of input arguments (should be 2 or 4)');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a Star');
            end
                          
            if lb >= 0
                S = Star(I.V, I.C, I.d, I.predicate_lb, I.predicate_ub);
            elseif ub <= 0
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                S = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub);
            elseif lb < 0 && ub > 0
                
                S1 = ub*(ub-lb)/2; % area of the first candidate abstract-domain
                S2 = -lb*(ub-lb)/2; % area of the second candidate abstract-domain  
                
                n = I.nVar + 1;
                                
                % constraint 1: y[index] = ReLU(x[index]) >= 0
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                
                
                % constraint 2: y[index] = ReLU(x[index]) >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                    
                % constraint 3: y[index] <= ub(x[index] -lb)/(ub - lb)
                C3 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                d3 = -(ub*lb/(ub-lb)) + ub*I.V(index, 1)/(ub-lb);
                               
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;
                
                if S1 < S2
                    % get first cadidate as resulted abstract-domain
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                    new_pred_lb = [I.predicate_lb; 0];
                    new_pred_ub = [I.predicate_ub; ub];
                    
                else
                    % choose the second candidate as the abstract-domain                                      
                    new_C = [C0; C2; C3];
                    new_d = [d0; d2; d3];
                    new_pred_lb = [I.predicate_lb; lb];
                    new_pred_ub = [I.predicate_ub; ub];
                                        
                end
                
                S = Star(new_V, new_C, new_d, new_pred_lb, new_pred_ub);
                
            end
                       
        end
        
        
        % over-approximate reachability analysis using abstract-domain
        function S = reach_abstract_domain(I)
            % @I: star input set
            % @S: star output set

            % author: Dung Tran
            % date: 4/3/2019


            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isempty(I)
                S = [];
            else    
                [lb, ub] = I.estimateRanges;
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    map = find(ub <= 0); % computation map
                    V = I.V;
                    V(map, :) = 0;
                    In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub);                    
                    map = find(lb <0 & ub > 0);
                    m = length(map); 
                    for i=1:m
                        fprintf('\nPerforming approximate PosLin_%d operation using Abstract Domain', map(i));
                        In = PosLin.stepReachAbstractDomain(In, map(i), lb(map(i)), ub(map(i)));
                    end
                    S = In;
             
                end
            
            end     
        
        end
        
    end
    
    
    
    methods(Static) % reachability analysis method using face-latice
        
        % future supporting method
        
    end
    
    
    methods(Static) % main reach method
        
        % main function for reachability analysis
        function R = reach(varargin)
            % @I: an array of star input sets
            % @method: 'exact-star' or 'approx-star' or 'approx-zono' or
            % 'abs-dom' or 'face-latice'
            % @option: = 'parallel' use parallel option
            %          = '' do use parallel option
            
            % author: Dung Tran
            % date: 27/2/2019
            
            switch nargin
                
                case 3
                    I = varargin{1};
                    method = varargin{2};
                    option = varargin{3};
                
                case 2
                    I = varargin{1};
                    method = varargin{2};
                    option = [];
                case 1
                    I = varargin{1};
                    method = 'exact-star';
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 1, 2 or 3)');
            end
            
            
            if strcmp(method, 'exact-star') % exact analysis using star
                
                R = PosLin.reach_star_exact_multipleInputs(I, option);
                
            elseif strcmp(method, 'exact-polyhedron') % exact analysis using polyhedron
                
                R = PosLin.reach_polyhedron_exact(I, option);
                
            elseif strcmp(method, 'approx-star')  % over-approximate analysis using star
                
                R = PosLin.reach_star_approx(I);
                
            elseif strcmp(method, 'approx-zono')  % over-approximate analysis using zonotope
                
                R = PosLin.reach_zono_approx(I);
                
            elseif strcmp(method, 'abs-dom')  % over-approximate analysis using abstract-domain
                
                R = PosLin.reach_abstract_domain(I);
                
            elseif strcmp(method, 'exact-face-latice') % exact analysis using face-latice
                fprintf('\nNNV have not yet support Exact Face-Latice Method');
            elseif strcmp(method, 'approx-face-latice') % over-approximate analysis using face-latice
                fprintf('\nNNV have not yet support Approximate Face-Latice Method');
            end
            
            
            
            
              
        end
        
        
    end
    
    
end

