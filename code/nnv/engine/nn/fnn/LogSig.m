classdef LogSig
    % LOGSIG Class contains methods for reachability analysis of layer with
    % Sigmoid activation function.
    % Reference: https://www.mathworks.com/help/deeplearning/ref/logsig.html
    % Author: Dung Tran
    % Date: 28/2/2019
    
    properties
        
    end
    
    methods(Static)  % evaluate method and over-approximate reachability analysis with stars
        
        % evaluation
        function y = evaluate(x)
            y = logsig(x);
        end
        
        
        % main method
        function S = reach_star_approx(varargin)
            % author: Dung Tran
            % date: 3/19/2020
            
            switch nargin
                case 1
                    I = varargin{1};
                    option = 'approx-star-no-split';
                case 2
                    I = varargin{1};
                    option = varargin{2};
                otherwise
                    error('Invalid number of input arguments');
            end
            
            if ~isa(I, 'Star')
                error('Input set is not a star set');
            end
            
            if strcmp(option, 'approx-star-no-split')
                S = LogSig.reach_star_approx_no_split(I);
            elseif strcmp(option, 'approx-star-split')
                S = LogSig.reach_star_approx_split(I);
            else
                error('Unknown reachability method');
            end
            
        end
        
        function S = reach_star_approx_multipleInputs(varargin)
            % author: Dung Tran
            % date: 3/27/2020
            
            switch nargin
                case 1
                    I = varargin{1};
                    parallel = []; % no parallel computation
                    method = 'approx-star-no-split';
                case 2
                    I = varargin{1};
                    parallel = varargin{2};
                    method = 'approx-star-no-split';
                case 3
                    I = varargin{1};
                    parallel = varargin{2};
                    method = varargin{3};
                otherwise
                    error('Invalid number of input arguments');
            end
            
            p = length(I);
            S = [];
            if isempty(parallel)
                
                for i=1:p
                    S =[S LogSig.reach_star_approx(I(i), method)];
                end
                
            elseif strcmp(parallel, 'parallel')
                
                parfor i=1:p
                    S =[S, LogSig.reach_star_approx(I(i), method)];
                end
                
            else
                error('Unknown parallel computation option');
            end

        end
        
        % reachability method with star
        function S = reach_star_approx_no_split(I)
            % @I: the input star set
            % @S: a star set output
            
            % author: Dung Tran
            % date: 3/19/2020
            
            B = I.getBox; 
            if isempty(B)
                S = [];
            else
                l = B.lb;
                u = B.ub;
                
                y_l = logsig(l);
                y_u = logsig(u);
                dy_l = logsig('dn', l);
                dy_u = logsig('dn', u);
                
                n = I.dim;
                S = I;
                for i=1:n
                    S = LogSig.stepLogSig_NoSplit(S, i, l(i), u(i), y_l(i), y_u(i), dy_l(i), dy_u(i)); 
                end
            end
            
            
        end
        
        % reachability method with star
        function S = reach_star_approx_split(I)
            % @I: the input star set
            % @S: an array of star set output
            
            % author: Dung Tran
            % date: 3/19/2020
            
            B = I.getBox; 
            if isempty(B)
                S = [];
            else
                l = B.lb;
                u = B.ub;
                
                y_l = logsig(l);
                y_u = logsig(u);
                dy_l = logsig('dn', l);
                dy_u = logsig('dn', u);
                
                n = I.dim;
                S = I;
                for i=1:n
                    m = length(S);
                    O = [];
                    for j=1:m
                        O = [O LogSig.stepLogSig_Split(S(j), i, l(i), u(i), y_l(i), y_u(i), dy_l(i), dy_u(i))];
                    end
                    S = O;
                end
            end
            
            
        end
        
        % stepLogSig
        function S = stepLogSig_Split(I, index, l, u, y_l, y_u, dy_l, dy_u)
            % @I: input star set
            % @index: index of the neuron
            % @l: l = min(x[index]), lower bound at neuron x[index] 
            % @u: u = min(x[index]), upper bound at neuron x[index]
            % @y_l: = logsig(l); output of logsig at lower bound
            % @y_u: = logsig(u); output of logsig at upper bound
            % @dy_l: derivative of LogSig at the lower bound
            % @dy_u: derivative of LogSig at the upper bound
            
            % @S: output star set
            
            % author: Dung Tran
            % date: 3/19/2020
            
            if l == u
                
                new_V = I.V;
                new_V(index, 1:I.nVar+1) = 0;
                new_V(index, 1) = y_l;
                S = Star(new_V, I.C, I.d, I.predicate_lb, I.predicate_ub);
                
            elseif l >= 0
                % y is convex when x >= 0
                % constraint 1: y <= y'(l) * (x - l) + y(l)
                % constarint 2: y <= y'(u) * (x - u) + y(u) 
                % constraint 3: y >= (y(u) - y(l)) * (x - l) / (u - l) + y(l);
                
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y <= y'(l) * (x - l) + y(l)
                C1 = [-dy_l*I.V(index, 2:n) 1];
                d1 = dy_l * I.V(index, 1) - dy_l*l + y_l; 
                % constraint 2: y <= y'(u) * (x - u) + y(u)
                C2 = [-dy_u*I.V(index, 2:n) 1];
                d2 = dy_u * I.V(index, 1) - dy_u*u + y_u;
                % constraint 3: y >= (y(u) - y(l)) * (x - l) / (u - l) + y(l);
                a = (y_u - y_l)/(u - l);
                C3 = [a*I.V(index, 2:n) -1];
                d3 = a*l - y_l - a*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1; 
                
                % update predicate bound
                new_predicate_lb = [I.predicate_lb; y_l]; 
                new_predicate_ub = [I.predicate_ub; y_u];
                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);
                
                
            elseif u <= 0
                % y is concave when x <= 0
                % constraint 1: y >= y'(l) * (x - l) + y(l)
                % constraint 2: y >= y'(u) * (x - u) + y(u)
                % constraint 3: y <= (y(u) - y(l)) * (x -l) / (u - l) + y(l);
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y >= y'(l) * (x - l) + y(l)
                C1 = [dy_l*I.V(index, 2:n) -1];
                d1 = -dy_l * I.V(index, 1) + dy_l*l - y_l; 
                % constraint 2: y >= y'(u) * (x - u) + y(u)
                C2 = [dy_u*I.V(index, 2:n) -1];
                d2 = -dy_u * I.V(index, 1) + dy_u*u - y_u;
                % constraint 3: y <= (y(u) - y(l)) * (x - l) / (u - l) + y(l);
                a = (y_u - y_l)/(u - l);
                C3 = [-a*I.V(index, 2:n) 1];
                d3 = -a*l + y_l + a*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1; 
                
                % update predicate bound
                new_predicate_lb = [I.predicate_lb; y_l]; 
                new_predicate_ub = [I.predicate_ub; y_u];
                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);
                
                
            elseif l <0 && u >0
                % y is concave for x in [l, 0] and convex for x
                % in [0, u]
                % split can be done here 
                
                % case 1: x in [l, 0]
                % y'(0) = 0.25
                % y is concave when x <= 0
                % constraint 1: y >= y'(l) * (x - l) + y(l)
                % constraint 2: y >= y'(0) * (x) + y(0)
                % constraint 3: y <= (y(0) - y(l)) * (x -l) / (0 - l) + y(l);
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y >= y'(l) * (x - l) + y(l)
                C1 = [dy_l*I.V(index, 2:n) -1];
                d1 = -dy_l * I.V(index, 1) + dy_l*l - y_l; 
                % constraint 2: y >= y'(0) * (x - 0) + y(0) = 0.25*x + 0.5
                C2 = [0.25*I.V(index, 2:n) -1];
                d2 = -0.25 * I.V(index, 1) - 0.5;
                % constraint 3: y <= (y(0) - y(l)) * (x - l) / (0 - l) + y(l);
                a = (0.5 - y_l)/(0 - l);
                C3 = [-a*I.V(index, 2:n) 1];
                d3 = -a*l + y_l + a*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1; 
                
                % update predicate bound
                new_predicate_lb = [I.predicate_lb; y_l]; 
                new_predicate_ub = [I.predicate_ub; 0.5];
                S1 = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);
                
                % case 2: x in [0, u] 
                % y is convex when x >= 0
                % constraint 1: y <= y'(0) * (x - 0) + y(0) = 0.25*x + 0.5
                % constarint 2: y <= y'(u) * (x - u) + y(u) 
                % constraint 3: y >= (y(u) - y(0)) * (x - 0) / (u - 0) + y(0);
                
                % over-approximation constraints 
                % constraint 1: y <= y'(0) * (x - 0) + y(0) = 0.25*x + 0.5
                C1 = [-0.25*I.V(index, 2:n) 1];
                d1 = 0.25 * I.V(index, 1) + 0.5; 
                % constraint 2: y <= y'(u) * (x - u) + y(u)
                C2 = [-dy_u*I.V(index, 2:n) 1];
                d2 = dy_u * I.V(index, 1) - dy_u*u + y_u;
                % constraint 3: y >= (y(u) - y(0)) * (x - 0) / (u - 0) + y(0);
                a = (y_u - 0.5)/u;
                C3 = [a*I.V(index, 2:n) -1];
                d3 = -0.5 - a*I.V(index, 1);
                
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1; 
                
                % update predicate bound
                new_predicate_lb = [I.predicate_lb; 0.5]; 
                new_predicate_ub = [I.predicate_ub; y_u];
                S2 = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);
                
                S = [S1 S2];
            end

        end
        
        % stepLogSig
        function S = stepLogSig_NoSplit(I, index, l, u, y_l, y_u, dy_l, dy_u)
            % @I: input star set
            % @index: index of the neuron
            % @l: l = min(x[index]), lower bound at neuron x[index] 
            % @u: u = min(x[index]), upper bound at neuron x[index]
            % @y_l: = logsig(l); output of logsig at lower bound
            % @y_u: = logsig(u); output of logsig at upper bound
            % @dy_l: derivative of LogSig at the lower bound
            % @dy_u: derivative of LogSig at the upper bound
            
            % @S: output star set
            
            % author: Dung Tran
            % date: 3/19/2020
            
            
            if l == u
                new_V = I.V;
                new_V(index, 1:I.nVar+1) = 0;
                new_V(index, 1) = y_l;
                S = Star(new_V, I.C, I.d, I.predicate_lb, I.predicate_ub);                
            elseif l >= 0
                % y is convex when x >= 0
                % constraint 1: y <= y'(l) * (x - l) + y(l)
                % constarint 2: y <= y'(u) * (x - u) + y(u) 
                % constraint 3: y >= (y(u) - y(l)) * (x - l) / (u - l) + y(l);
                
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y <= y'(l) * (x - l) + y(l)
                C1 = [-dy_l*I.V(index, 2:n) 1];
                d1 = dy_l * I.V(index, 1) - dy_l*l + y_l; 
                % constraint 2: y <= y'(u) * (x - u) + y(u)
                C2 = [-dy_u*I.V(index, 2:n) 1];
                d2 = dy_u * I.V(index, 1) - dy_u*u + y_u;
                % constraint 3: y >= (y(u) - y(l)) * (x - l) / (u - l) + y(l);
                a = (y_u - y_l)/(u - l);
                C3 = [a*I.V(index, 2:n) -1];
                d3 = a*l - y_l - a*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1; 
                
                % update predicate bound
                new_predicate_lb = [I.predicate_lb; y_l]; 
                new_predicate_ub = [I.predicate_ub; y_u];
                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);
                
                
            elseif u <= 0
                % y is concave when x <= 0
                % constraint 1: y >= y'(l) * (x - l) + y(l)
                % constraint 2: y >= y'(u) * (x - u) + y(u)
                % constraint 3: y <= (y(u) - y(l)) * (x -l) / (u - l) + y(l);
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y >= y'(l) * (x - l) + y(l)
                C1 = [dy_l*I.V(index, 2:n) -1];
                d1 = -dy_l * I.V(index, 1) + dy_l*l - y_l; 
                % constraint 2: y >= y'(u) * (x - u) + y(u)
                C2 = [dy_u*I.V(index, 2:n) -1];
                d2 = -dy_u * I.V(index, 1) + dy_u*u - y_u;
                % constraint 3: y <= (y(u) - y(l)) * (x - l) / (u - l) + y(l);
                a = (y_u - y_l)/(u - l);
                C3 = [-a*I.V(index, 2:n) 1];
                d3 = -a*l + y_l + a*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1; 
                
                % update predicate bound
                new_predicate_lb = [I.predicate_lb; y_l]; 
                new_predicate_ub = [I.predicate_ub; y_u];
                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);
                
                
            elseif l <0 && u >0
                % y is concave for x in [l, 0] and convex for x
                % in [0, u]
                % split can be done here 
                
                % combine two stars into one star               
                x1 = (0.5 + dy_u * u - y_u)/(dy_u - 0.25);
                y1 = 0.25 * x1 + 0.5;
                x2 = (0.5 + dy_l * l - y_l)/(dy_l - 0.25);
                y2 = 0.25 * x2 + 0.5; 
                
                % over-approximation constraints 
                % constraint 1: y >= y'(l) * (x - l) + y(l)
                % constraint 2: y <= y'(u) * (x - u) + y(u)
                % constraint 3: y <= (y(x1) - y(l))*(x - l)/(x1 - l) + y(l);
                % constraint 4: y >= (y(x2) - y(u)) * (x - u)/(x2 - u) + y(u)
                
                               
                n = I.nVar + 1;
                
                % constraint 1: y >= y'(l) * (x - l) + y(l)
                C1 = [dy_l*I.V(index, 2:n) -1];
                d1 = -dy_l * I.V(index, 1) + dy_l*l - y_l; 
                % constraint 2: y <= y'(u) * (x - u) + y(u)
                C2 = [-dy_u*I.V(index, 2:n) 1];
                d2 = dy_u * I.V(index, 1) - dy_u*u + y_u;
                % constraint 3: y <= (y(x1) - y(l))*(x - l)/(x1 - l) + y(l);
                a = (y1 - y_l)/(x1 - l);
                C3 = [-a*I.V(index, 2:n) 1];
                d3 = -a*l + y_l + a*I.V(index, 1);
                                
                % constraint 4: y >= (y(x2) - y(u)) * (x - u)/(x2 - u) + y(u);
                a = (y2 - y_u)/(x2 - u);
                C4 = [a*I.V(index, 2:n) -1];
                d4 = a*u - y_u - a*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3; C4];
                new_d = [d0; d1; d2; d3; d4];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1; 
                
                % update predicate bound
                new_predicate_lb = [I.predicate_lb; y_l]; 
                new_predicate_ub = [I.predicate_ub; y_u];
                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);
                                
            end

        end
          
    end
    
    
    methods(Static) % over-approximate reachability analysis using Zonotope
        
        function Z = reach_zono_approx(I)
            % @I: zonotope input set
            % @Z: zonotope output set
            
            % author: Dung Tran
            % date: 5/3/2019
            
            % reference: Fast and Effective Robustness Certification,
            % Gagandeep Singh, NIPS, 2018
            
            if ~isa(I, 'Zono')
                error('Input set is not a Zonotope');
            end
            
            B = I.getBox;
            
            lb = B.lb;
            ub = B.ub;
            G = [logsig('dn', lb) logsig('dn', ub)];
            gamma_opt = min(G, [], 2);
            gamma_mat = diag(gamma_opt);
            mu1 = 0.5 * (logsig(ub) + logsig(lb) - gamma_mat * (ub + lb));
            mu2 = 0.5 * (logsig(ub) - logsig(lb) - gamma_mat * (ub - lb));
            Z1 = I.affineMap(gamma_mat, mu1);
            new_V = diag(mu2);
            V = [Z1.V new_V];
            Z = Zono(Z1.c, V);
            
        end
        
        % dealing with multiple inputs in parallel
        function S = reach_zono_approx_multipleInputs(varargin)
            % author: Dung Tran
            % date: 3/27/2020
            
            switch nargin
                case 1
                    I = varargin{1};
                    parallel = []; % no parallel computation
                case 2
                    I = varargin{1};
                    parallel = varargin{2};
                otherwise
                    error('Invalid number of input arguments, should be 1 or 2');
            end
            
            p = length(I);
            S = [];
            if isempty(parallel)
                
                for i=1:p
                    S =[S LogSig.reach_zono_approx(I(i))];
                end
                
            elseif strcmp(parallel, 'parallel')
                
                parfor i=1:p
                    S =[S, LogSig.reach_zono_approx(I(i))];
                end
                
            else
                error('Unknown parallel computation option');
            end

        end
        
    end
    
methods(Static) % over-approximate reachability analysis using abstract domain
    
    
    function S = stepLogSig_absdom(I, index, l, u, y_l, y_u, dy_l, dy_u)
        % @I: input star set
        % @index: index of the neuron
        % @l: l = min(x[index]), lower bound at neuron x[index] 
        % @u: u = min(x[index]), upper bound at neuron x[index]
        % @y_l: = logsig(l); output of logsig at lower bound
        % @y_u: = logsig(u); output of logsig at upper bound
        % @dy_l: derivative of LogSig at the lower bound
        % @dy_u: derivative of LogSig at the upper bound

        % @S: output star set

        % author: Dung Tran
        % date: 3/27/2020
                
        if l == u
            new_V = I.V;
            new_V(index, 1:I.nVar+1) = 0;
            new_V(index, 1) = y_l;
            S = Star(new_V, I.C, I.d, I.predicate_lb, I.predicate_ub);                
        elseif l >= 0
            % y is convex when x >= 0
            % constraint 2: y <= y'(u) * (x - u) + y(u) 
            % constraint 3: y >= (y(u) - y(l)) * (x - l) / (u - l) + y(l);


            n = I.nVar + 1;
            % over-approximation constraints 
            % constraint 2: y <= y'(u) * (x - u) + y(u)
            C2 = [-dy_u*I.V(index, 2:n) 1];
            d2 = dy_u * I.V(index, 1) - dy_u*u + y_u;
            % constraint 3: y >= (y(u) - y(l)) * (x - l) / (u - l) + y(l);
            a = (y_u - y_l)/(u - l);
            C3 = [a*I.V(index, 2:n) -1];
            d3 = a*l - y_l - a*I.V(index, 1);

            m = size(I.C, 1);
            C0 = [I.C zeros(m, 1)];
            d0 = I.d;
            new_C = [C0; C2; C3];
            new_d = [d0; d2; d3];
            new_V = [I.V zeros(I.dim, 1)];
            new_V(index, :) = zeros(1, n+1);
            new_V(index, n+1) = 1; 

            % update predicate bound
            new_predicate_lb = [I.predicate_lb; y_l]; 
            new_predicate_ub = [I.predicate_ub; y_u];
            S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);


        elseif u <= 0
            % y is concave when x <= 0
            % constraint 1: y >= y'(l) * (x - l) + y(l)
            % constraint 3: y <= (y(u) - y(l)) * (x -l) / (u - l) + y(l);

            n = I.nVar + 1;
            % over-approximation constraints 
            % constraint 1: y >= y'(l) * (x - l) + y(l)
            C1 = [dy_l*I.V(index, 2:n) -1];
            d1 = -dy_l * I.V(index, 1) + dy_l*l - y_l; 
            % constraint 3: y <= (y(u) - y(l)) * (x - l) / (u - l) + y(l);
            a = (y_u - y_l)/(u - l);
            C3 = [-a*I.V(index, 2:n) 1];
            d3 = -a*l + y_l + a*I.V(index, 1);

            m = size(I.C, 1);
            C0 = [I.C zeros(m, 1)];
            d0 = I.d;
            new_C = [C0; C1; C3];
            new_d = [d0; d1; d3];
            new_V = [I.V zeros(I.dim, 1)];
            new_V(index, :) = zeros(1, n+1);
            new_V(index, n+1) = 1; 

            % update predicate bound
            new_predicate_lb = [I.predicate_lb; y_l]; 
            new_predicate_ub = [I.predicate_ub; y_u];
            S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);


        elseif l <0 && u >0
            % y is concave for x in [l, 0] and convex for x
            % in [0, u]
            % split can be done here 

            % combine two stars into one star               
            x1 = (0.5 + dy_u * u - y_u)/(dy_u - 0.25);
            y1 = 0.25 * x1 + 0.5;
            x2 = (0.5 + dy_l * l - y_l)/(dy_l - 0.25);
            y2 = 0.25 * x2 + 0.5; 

            % over-approximation constraints 
            % constraint 1: y >= y'(l) * (x - l) + y(l)
            % constraint 2: y <= y'(u) * (x - u) + y(u)
           
            n = I.nVar + 1;

            dy_min = min(dy_l, dy_u);
            % constraint 1: y >= y'_min * (x - l) + y(l)
            C1 = [dy_min*I.V(index, 2:n) -1];
            d1 = -dy_min * I.V(index, 1) + dy_min*l - y_l; 
            % constraint 2: y <= y'_min * (x - u) + y(u)
            C2 = [-dy_min*I.V(index, 2:n) 1];
            d2 = dy_min * I.V(index, 1) - dy_min*u + y_u;
            
            m = size(I.C, 1);
            C0 = [I.C zeros(m, 1)];
            d0 = I.d;
            new_C = [C0; C1; C2];
            new_d = [d0; d1; d2];
            new_V = [I.V zeros(I.dim, 1)];
            new_V(index, :) = zeros(1, n+1);
            new_V(index, n+1) = 1; 

            % update predicate bound
            new_predicate_lb = [I.predicate_lb; y_l]; 
            new_predicate_ub = [I.predicate_ub; y_u];
            S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);

        end

        
        
    end
    
    function S = reach_absdom_approx(I)
        % @I: star input set
        % @Z: Star output set

        % author: Dung Tran
        % date: 3/27/2020

        % reference: An abstract domain for certifying neural networks. Proceedings of the ACM on Programming Languages,
        % Gagandeep Singh, POPL, 2019

        if ~isa(I, 'Star')
            error('Input set is not a Star');
        end

        [l, u] = I.estimateRanges;  

        y_l = logsig(l);
        y_u = logsig(u);
        dy_l = logsig('dn', l);
        dy_u = logsig('dn', u);

        n = I.dim;
        S = I;
        for i=1:n
            S = LogSig.stepLogSig_absdom(S, i, l(i), u(i), y_l(i), y_u(i), dy_l(i), dy_u(i)); 
        end

    end
    
    % dealing with multiple inputs in parallel
    function S = reach_absdom_approx_multipleInputs(varargin)
        % author: Dung Tran
        % date: 3/27/2020

        switch nargin
            case 1
                I = varargin{1};
                parallel = []; % no parallel computation
            case 2
                I = varargin{1};
                parallel = varargin{2};
            otherwise
                error('Invalid number of input arguments, should be 1 or 2');
        end

        p = length(I);
        S = [];
        if isempty(parallel)

            for i=1:p
                S =[S LogSig.reach_absdom_approx(I(i))];
            end

        elseif strcmp(parallel, 'parallel')

            parfor i=1:p
                S =[S, LogSig.reach_absdom_approx(I(i))];
            end

        else
            error('Unknown parallel computation option');
        end

    end

end


methods(Static) % main reach method
    
    % main function for reachability analysis
    function R = reach(varargin)
        % @I: an array of star input sets
        % @method: 'approx-star' or 'approx-zono' or 'abs-dom' 
        % @option: = 'parallel' or [] using parallel computation or not

        % author: Dung Tran
        % date: 3/27/2019

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
                method = 'approx-star';
                option = [];
            otherwise
                error('Invalid number of input arguments (should be 1, 2 or 3)');
        end


        if strcmp(method, 'approx-star') % exact analysis using star

            R = LogSig.reach_star_approx_multipleInputs(I, option);

        elseif strcmp(method, 'approx-zono')  % over-approximate analysis using zonotope

            R = LogSig.reach_zono_approx_multipleInputs(I, option);

        elseif strcmp(method, 'abs-dom')  % over-approximate analysis using abstract-domain

            R = LogSig.reach_absdom_approx_multipleInputs(I, option);

        else
            error('Unknown or unsupported reachability method for layer with LogSig activation function');
        end

    end
    
end



    
end

