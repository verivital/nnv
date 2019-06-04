classdef SatLin
    % SATLIN : class for computing reachable set of Satlin Transfer Function 
    %   Reference: https://www.mathworks.com/help/deeplearning/ref/satlin.html
    % Author: Dung Tran
    % Date: 27/2/2019
    
    properties
        
    end
    
    methods(Static) % evaluate method and reachability analysis with stars    
        
        % evaluation
        function y = evaluate(x)
            y = satlin(x);
        end
            
        % stepReach method, compute reachable set for a single step
        function S = stepReach(I, index)
            % @I: single star set input
            % @index: index of the neural performing stepSatLin
            % @S: star output set
            
            % author: Dung Tran
            % date: 27/2/2019
            
            
            if ~isa(I, 'Star')
                error('Input is not a star set');
            end
            
            % case 1: x(index) <= 0, SatLin(x(index)) = 0
            C0 = I.V(index, 2:I.nVar + 1);
            d0 = I.V(index, 1);
            C1 = [I.C; C0];
            d1 = [I.d; -d0];
            V1 = I.V;
            V1(index, :) = zeros(1, I.nVar + 1);
            S1 = Star(V1, C1, d1);
            if S1.isEmptySet
                S1 = [];
            end
            
            % case 2: 0 <= x(index) <= 1, SatLin(x(index)) = x(index)            
            V2 = I.V;
            C2 = [I.C; C0; -C0];
            d2 = [I.d; 1-d0; d0];
            S2 = Star(V2, C2, d2);
            if S2.isEmptySet
                S2 = [];
            end
            
            % case 3: x(index) >= 1
            C3 = [I.C; -C0];
            d3 = [I.d; d0 - 1];
            V3 = I.V;
            V3(index, 1) = 1;
            V3(index, 2:I.nVar + 1) = zeros(1, I.nVar);
            S3 = Star(V3, C3, d3);
            if S3.isEmptySet
                S3 = [];
            end
            
            S = [S1 S2 S3];
                     
        end
        
        
        % stepReach with multiple inputs
        function S = stepReachMultipleInputs(varargin)
            % @I: an array of stars
            % @index: index where stepReach is performed
            % @option: = 'parallel' use parallel computing
            %          = not declare -> don't use parallel computing
            
            % author: Dung Tran
            % date: 27/2/2019
            
            switch nargin
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    option = varargin{3};
                case 2
                    I = varargin{1};
                    index = varargin{2};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2 or 3)');
            end
            
            
            
            p = length(I);
            S = [];
            
            if isempty(option)
                
                for i=1:p
                    S =[S, SatLin.stepReach(I(i), index)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    S =[S, SatLin.stepReach(I(i), index)];
                end
                
            else
                error('Unknown option');
            end
            
            
        end
       
        
        
        % exact reachability analysis using Star
        function S = reach_star_exact(I, option)
            % @I: an array of star input sets
            % @option: = 'parallel' use parallel option
            %          = '' do use parallel option
            
            % author: Dung Tran
            % date: 27/2/2019
            
            if ~isempty(I)       
                dim = I(1).dim;
                In = I;
                for i=1:dim
                    fprintf('\nPerforming exact SatLin_%d operation using Star', i);
                    In = SatLin.stepReachMultipleInputs(In, i, option);
                end             
                
                S = In;
            else
                S = [];
            end
            
              
        end
        
        
        % step over-approximate reachability analysis using Star
        function S = stepReachStarApprox(I, index)
            % @I: star input set
            % @index: index of the neuron where we perform step Reach
            % @S: Star output set 
            
            % author: Dung Tran
            % date: 4/3/2019
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end
            
            [lb, ub] = I.getRange(index);
            
            if ub <= 0
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                S = Star(V, I.C, I.d);
            end
            if lb >= 1
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                V(index, 1) = 1;
                S = Star(V, I.C, I.d);                
            end
            if (1 > lb) && (lb > 0) && ub > 1
                % constraint 1: y(index) <= x[index]
                C1 = [-I.V(index, 2:I.nVar + 1) 1];
                d1 = I.V(index, 1);
                % constraint 2: y[index] <= 1
                C2 = zeros(1, I.nVar + 1);
                C2(I.nVar + 1) = 1;
                d2 = 1;
                % constraint 3: y[index] >= (1-lb)x/(ub-lb) + lb*(ub-1)/(ub-lb)
                C3 = [((1-lb)/(ub-lb))*I.V(index, 2:I.nVar+1) -1];
                d3 = -lb*(ub-1)/(ub-lb) - (1-lb)*I.V(index,1)/(ub-lb);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, I.nVar+2);
                new_V(index, I.nVar+2) = 1;

                S = Star(new_V, new_C, new_d);
            end
            
            if lb >= 0 && ub <= 1
                S = Star(I.V, I.C, I.d);
            end
            if lb < 0 && (0 < ub) && (ub <= 1)
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y[index] = >= 0
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                % constraint 2: y[index] >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                % constraint 3: y[index] <= ub(x[index] -lb)/(ub - lb)
                C3 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                d3 = -(ub*lb/(ub-lb))*(1 - I.V(index, 1));

                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;

                S = Star(new_V, new_C, new_d);               
            end
            
            if lb < 0 && ub > 1
                % constraint 1: y[index] >= 0
                n = I.nVar + 1;
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                % constraint 2: y[index] <= 1
                C2 = zeros(1, n);
                C2(n) = 1;
                d2 = 1;
                % constraint 3: y[index] <= x/(1 -lb) - lb/(1-lb)
                C3 = [(-1/(1-lb))*I.V(index, 2:n) 1];
                d3 = (1/(1-lb))*I.V(index, 1) - lb/(1-lb);
                % constraint 4: y[index] >=  x/ub
                C4 = [(1/ub)*I.V(index, 2:n) -1];
                d4 = -(1/ub)*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3; C4];
                new_d = [d0; d1; d2; d3; d4];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;

                S = Star(new_V, new_C, new_d); 
                
            end
                 
        end
        
        % over-approximate reachability analysis using star
        function S = reach_star_approx(I)
            % @I: star input set
            % @S: star output set
            
            % author: Dung Tran
            % date: 4/3/2019
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isEmptySet(I)
                S = [];
            else
                In = I;
                for i=1:I.dim
                    fprintf('\nPerforming exact SatLin_%d operation using star', i);
                    In = SatLin.stepReachStarApprox(In, i);
                end
                S = In;
            end
        end
        
             
    end
    
    % exact analysis using polyhedron
    methods(Static)
        
        % step reachability analysis using polyhedron
        function P = stepReachPolyhedronExact(I, index)
            % @I: polyhedron input set
            % @index: index of the neuron performing step reach
            % @P: polyhedron output set 
            
            % author: Dung Tran
            % date: 6/4/2019
            
            
            if ~isa(I, 'Polyhedron')
                error('Input is not a polyhedron');
            end
            
            dim = I.Dim;
            
            
            
            Im = eye(dim);
            Im(index, index) = 0;
            % case 1: x(index) <= 0, SatLin(x(index)) = 0
            Z = zeros(1, dim);
            Z(1, index) = 1;
            A = vertcat(I.A, Z);
            b = [I.b; 0];
            R1 = Polyhedron('A', A, 'b', b, 'Ae', I.Ae, 'be', I.be);
            R1 = R1.affineMap(Im, 'vrep');
            
            if R1.isEmptySet
                R1 = [];
            end
            
            % case 2: x[index] > 0 , x[index] < 1 -> SatLin(x[index]) = x[index]          
            Z = [zeros(1, dim); zeros(1, dim)];
            Z(1, index) = -1;
            Z(2, index) = 1; 
            A = vertcat(I.A, Z);
            b = [I.b; 0; 1];
            R2 = Polyhedron('A', A, 'b', b, 'Ae', I.Ae, 'be', I.be);
            
            if R2.isEmptySet
               R2 = [];
            end

            % case 3: x[index] >= 1
            Z = zeros(1, dim);
            Z(1, index) = -1;
            A = vertcat(I.A, Z);
            b = [I.b; -1];
            R3 = Polyhedron('A', A, 'b', b, 'Ae', I.Ae, 'be', I.be);
            R3 = R3.affineMap(Im, 'vrep');
            z = zeros(dim, 1);
            z(index) = 1;
            
            if ~R3.isEmptySet
                R3 = R3 + z;
            else
                R3 = [];
            end
            
            P = [R1 R2 R3];
            
                              
        end
        
        % stepReachPolyhedron with multiple inputs
        function P = stepReachPolyhedronMultipleInputs(varargin)
            % @I: an array of stars
            % @index: index where stepReach is performed
            % @option: = 'parallel' use parallel computing
            %          = not declare -> don't use parallel computing
            
            % author: Dung Tran
            % date: 6/4/2019
            
            switch nargin
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    option = varargin{3};
                case 2
                    I = varargin{1};
                    index = varargin{2};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2 or 3)');
            end
            
            
            
            p = length(I);
            P = [];
            
            if isempty(option)
                
                for i=1:p
                    P =[P, SatLin.stepReachPolyhedronExact(I(i), index)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    P =[P, SatLin.stepReachPolyhedronExact(I(i), index)];
                end
                
            else
                error('Unknown option');
            end
            
            
        end
        
        
        % exact reachability analysis using polyhedorn
        function S = reach_polyhedron_exact(I, option)
            % @I: an array of star input sets
            % @option: = 'parallel' use parallel option
            %          = '' do use parallel option
            
            % author: Dung Tran
            % date: 6/4/2019
            
            if ~isempty(I)       
                dim = I(1).Dim;
                In = I;
                for i=1:dim
                    fprintf('\nPerforming exact SatLin_%d operation using Polyhedron', i);
                    In = SatLin.stepReachPolyhedronMultipleInputs(In, i, option);
                end             
                
                S = In;
            else
                S = [];
            end
            
              
        end
           
    end

    
    methods(Static) % over-approximate reachability analysis using zonotope
        
        % step reachability analysis using zonotope
        function Z = stepReachZonoApprox(I, index)
            % @I: zonotope input set
            % @index: index of the neuron performing step reach
            % @Z: zonotope output set
            
            % author: Dung Tran
            % date: 5/3/2019
            
            if ~isa(I, 'Zono')
                error('Input set is not a Zonotope');
            end
            
            [lb, ub] = I.getRange(index);
              
            if ub <= 0
                V = I.V;
                V(index, :) = zeros(1, size(V, 2));
                c = I.c;
                c(index) = 0;
                Z = Zono(c, V);
            end
            if lb >= 1
                V = I.V;
                V(index, :) = zeros(1, size(V,2));
                c = I.c;
                c(index) = 1;
                Z = Zono(c, V);                
            end
            
            if lb >= 0 && ub <= 1
                Z = Zono(I.c, I.V);
            end
            
            if (1 > lb) && (lb > 0) && ub > 1
                c = I.c;
                c(index) = c(index) + (1-ub)/2;
                V = zeros(I.dim, 1);
                V(index) = (1-ub)/2;
                V = [I.V V];
                Z = Zono(c, V);
            end
            
            if lb < 0 && (0 < ub) && (ub <= 1)
                
                c = I.c;
                lamda_opt = ub/(ub - lb);
                mu = -0.5*ub*lb/(ub - lb);
                c(index) = lamda_opt * c(index) + mu;
                V = zeros(I.dim, 1);
                V(index) = mu;
                V = [I.V V];
                Z = Zono(c, V);
                     
            end
            
            if lb < 0 && ub > 1
                
                % x + 1 -ub <= y <= x - lb, lb <= x <= ub
                c = I.c;
                mu = (1 + lb - ub)/2;
                c(index) = c(index) - lb  + mu;
                V = zeros(I.dim, 1);
                V(index) = mu;
                V = [I.V V];
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
            for i=1:I.dim
                fprintf('\nPerforming approximate SatLin_%d operation using Zonotope', i);
                In = SatLin.stepReachZonoApprox(In, i);
            end
            Z = In;
            
        end
        
        
    end
    
    
    methods(Static) % reachability analysis using abstract-domain
        
        % step over-approximate reachability analysis using abstract-domain
        % we extend abstract-domain method from ReLU to SatLin
        % we use star set to represent abstract-domain
        function S = stepReachAbstractDomain(I, index)
            % @I: star-input set
            % @index: index of neuron performing stepReach
            % @S: star output set represent abstract-domain of the output
            % set
            
            % author: Dung Tran
            % date: 16/3/2019
        
            % reference: An Abstract Domain for Certifying Neural Networks,
            % Gagandeep Singh, POPL 2019
           
            if ~isa(I, 'Star')
                error('Input is not a Star');
            end
               
            [lb, ub] = I.getRange(index);
            
            if ub <= 0
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                S = Star(V, I.C, I.d);
            end
            if lb >= 1
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                V(index, 1) = 1;
                S = Star(V, I.C, I.d);                
            end
            if (1 > lb) && (lb > 0) && ub > 1
                % constraint 1: y(index) <= x[index]
                C1 = [-I.V(index, 2:I.nVar + 1) 1];
                d1 = I.V(index, 1);
                % constraint 2: y[index] <= 1
                C2 = zeros(1, I.nVar + 1);
                C2(I.nVar + 1) = 1;
                d2 = 1;
                % constraint 3: y[index] >= (1-lb)x/(ub-lb) + lb*(ub-1)/(ub-lb)
                C3 = [((1-lb)/(ub-lb))*I.V(index, 2:I.nVar+1) -1];
                d3 = -lb*(ub-1)/(ub-lb) - (1-lb)*I.V(index,1)/(ub-lb);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, I.nVar+2);
                new_V(index, I.nVar+2) = 1;
                
                S1 = (1-lb)*(ub-lb)/2; % area of the first candidate abstract-domain
                S2 = (ub-lb)*(ub-1)/2;% area of the second candidate abstract-domain
                
                if S1 < S2
                    % get first cadidate as resulted abstract-domain
                    new_C = [C0; C2; C3];
                    new_d = [d0; d2; d3];
                else
                    % get second candidate as resulted abstract-domain
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                end

                S = Star(new_V, new_C, new_d);
            end
            
            if lb >= 0 && ub <= 1
                S = Star(I.V, I.C, I.d);
            end
            if lb < 0 && (0 < ub) && (ub <= 1)
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y[index] = >= 0
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                % constraint 2: y[index] >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                % constraint 3: y[index] <= ub(x[index] -lb)/(ub - lb)
                C3 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                d3 = -(ub*lb/(ub-lb))*(1 - I.V(index, 1));

                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;
                
                
                S1 = ub*(ub-lb)/2; % area of the first candidate abstract-domain
                S2 = -lb*(ub-lb)/2; % area of the second candidate abstract-domain
                
                if S1 < S2
                    % get first cadidate as resulted abstract-domain
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                else
                    % get second candidate as resulted abstract-domain
                    new_C = [C0; C2; C3];
                    new_d = [d0; d2; d3];
                end

                S = Star(new_V, new_C, new_d);               
            end
            
            if lb < 0 && ub > 1
                % constraint 1: y[index] >= 0
                n = I.nVar + 1;
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                % constraint 2: y[index] <= 1
                C2 = zeros(1, n);
                C2(n) = 1;
                d2 = 1;
                % constraint 3: y[index] <= x/(1 -lb) - lb/(1-lb)
                C3 = [(-1/(1-lb))*I.V(index, 2:n) 1];
                d3 = (1/(1-lb))*I.V(index, 1) - lb/(1-lb);
                % constraint 4: y[index] >=  x/ub
                C4 = [(1/ub)*I.V(index, 2:n) -1];
                d4 = -(1/ub)*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;
                
                S1 = (ub-lb)*(ub-lb)/(2*(1-lb)); % area of the first candidate abstract domain
                S2 = (1 - lb/ub)*(ub-lb)/2; % area of the second candidate abstract domain
                
                if S1 < S2
                    % get first cadidate as resulted abstract-domain
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                else
                    % get second candidate as resulted abstract-domain
                    new_C = [C0; C2; C4];
                    new_d = [d0; d2; d4];
                end
                
                S = Star(new_V, new_C, new_d); 
                
            end
                       
        end
        
        
        % over-approximate reachability analysis using abstract-domain
        function S = reach_abstract_domain(I)
            % @I: star input set
            % @S: star output set

            % author: Dung Tran
            % date: 16/3/2019


            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isEmptySet(I)
                S = [];
            else
                In = I;
                for i=1:I.dim
                    fprintf('\nPerforming approximate SatLin_%d operation using Abstract-Domain', i);
                    In = SatLin.stepReachAbstractDomain(In, i);
                end
                S = In;
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
                
                R = SatLin.reach_star_exact(I, option);
                
            elseif strcmp(method, 'exact-polyhedron') % exact analysis using star
                
                R = SatLin.reach_polyhedron_exact(I, option);
                
            elseif strcmp(method, 'approx-star')  % over-approximate analysis using star
                
                R = SatLin.reach_star_approx(I);
                
            elseif strcmp(method, 'approx-zono')  % over-approximate analysis using zonotope
                
                R = SatLin.reach_zono_approx(I);
                
            elseif strcmp(method, 'abs-dom')  % over-approximate analysis using abstract-domain
                
                R = SatLin.reach_abstract_domain(I);
                
            elseif strcmp(method, 'exact-face-latice') % exact analysis using face-latice
                fprintf('\nNNV have not yet support Exact Face-Latice Method');
            elseif strcmp(method, 'approx-face-latice') % over-approximate analysis using face-latice
                fprintf('\nNNV have not yet support Approximate Face-Latice Method');
            end
                        
        end
        
    end
        
    
end

