classdef RStar0
    % Relaxed Star zero constrained class
    % Sung Woo Choi: 12/17/2020
    
    % reference: An Abstract Domain for Certifying Neural Networks,
    % Gagandeep Singh, POPL 2019
    
    properties
        lower_a = {[]}; % a set of matrix for lower polyhedral constraint
        upper_a = {[]}; % a set of matrix for upper polyhedral constraint
        lb = {[]}; % a set of matrix for lower bound
        ub = {[]}; % a set of matrix for upper bound

        iter = inf; % number of iterations for back substitution
        dim = 0; % dimension of current relaxed polyhedron
    end
    
    methods
        
        % constructor
        function obj = RStar0(varargin)
            % @lower_a: a set of matrix for lower polyhedral constraint
            % @upper_a: a set of matrix for upper polyhedral constraint

            switch nargin
                
                case 5
                    lower_a = varargin{1};
                    upper_a = varargin{2};
                    lb = varargin{3};
                    ub = varargin{4};
                    iter = varargin{5};
                    
%                     [nL, mL] = size(lower_a);
%                     [nU, mU] = size(upper_a);
%                     [nl, ml] = size(lb);
%                     [nu, mu] = size(ub);
%                     
%                     if nL ~= nU || mL ~= mU
%                         error('Inconsistency between upper and lower polyhedral constraaints');
%                     end
%                     
%                     if ml ~= 1 || mu ~= 1
%                        error(' 
%                     end

%                     [n_nVar, m_nVar] = size(nVar);
%                     [nL, mL] = size(lower_a);
%                     [nU, mU] = size(upper_a);
%                     
%                     if mL ~= mU || nL ~= nU
%                         error('sizes of matrix of upper and lower polyhedral constraints do not match');
%                     elseif m_nVar ~= 1
%                         error('nVar must be one column vector');   
%                     end
                    
                    obj.lower_a = lower_a;
                    obj.upper_a = upper_a;
                    obj.lb = lb;
                    obj.ub = ub;
                    
                    obj.iter = iter;
                    len = length(lower_a);
                    obj.dim = size(lower_a{len}, 1);
                    
                case 4
                    lower_a = varargin{1};
                    upper_a = varargin{2};
                    lb = varargin{3};
                    ub = varargin{4};
                    iter = inf;
                    
                    obj = RStar0(lower_a, upper_a, lb, ub, iter);

                case 3
                    % construct RStar0 from lower and upper bound vector
                    lb = varargin{1};
                    ub = varargin{2};
                    iter = varargin{3};
                    
                    dim = size(lb, 1);
                    obj.lower_a{1} = [zeros(dim, 1) -eye(dim)]; 
                    obj.upper_a{1} = [zeros(dim, 1) eye(dim)];
                    obj.lb{1} = lb;
                    obj.ub{1} = ub;
                    
                    obj.iter = iter;
                    obj.dim = dim;
                case 2
                    I = varargin{1};
                    iter = varargin{2};
                    if isa(I, 'Polyhedron')
                        dim = I.Dim;
                        I.outerApprox;
                        l = I.Internal.lb;
                        u = I.Internal.ub;
                    elseif isa(I, 'Star')
                        dim = I.dim;
                        if ~isempty(I.state_lb) && ~isempty(I.state_ub)
                            l = I.state_lb;
                            u = I.state_ub;
                        else
                            [l, u] = I.getRanges();
                        end
                    elseif isa(I, 'Zono')
                        dim = I.dim;
                        [l, u] = I.getRanges();
                    else
                        error('Unkown imput set');
                    end
                    
                    if iter <= 0
                        error('Iteration must be greater than zero');
                    end

                    lower_a{1} = [zeros(dim, 1) eye(dim)];
                    upper_a{1} = [zeros(dim, 1) eye(dim)];
                    lb{1} = l;
                    ub{1} = u;
                    obj = RStar0(lower_a, upper_a, lb, ub, iter);
                
                case 1
                    I = varargin{1};
                    n = length(I);
                    if n > 1
                        for i = 1:n
                            obj(i) = RStar0(I(i), inf);
                        end
                    else
                        obj = RStar0(I, inf);
                    end
                    
                case 0
                    % create empty RStar0 seet (for preallocaation an array
                    % of RStar0)
                    obj.lower_a = {[]};
                    obj.upper_a = {[]};
                    obj.lb = {[]};
                    obj.ub = {[]};
                    obj.iter = inf;
                    obj.dim = 0;
                    
                otherwise
                    error('Invalid number of input arguments (should be 0 or 2 or 5)');
            end
        end
        
        % affine abstract mapping of RStar0 set
        function R = affineMap(varargin)
            switch nargin
                case 3
                    obj = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                case 2
                    obj = varargin{1};
                    W = varargin{2};
                    b = zeros(size(W,1),1);
            end
            
            [nW, mW] = size(W);
            [nb, mb] = size(b);
            
            if mW ~= obj.dim
                error('Inconsistency between the affine mapping matrix and dimension of the RStar0 set');
            end
            
            if mb > 1
                error('bias vector must be one column');
            end
            
            if nW ~= nb && nb~=0
                error('Inconsistency between the affine mapping matrix and the bias vector');
            end
            
            if nb == 0
                b = zeros(nW, 1);
            end

            % new lower and upper polyhedral contraints
            lower_a = obj.lower_a;
            upper_a = obj.upper_a;
            
            len = length(lower_a);
            lower_a{len+1} = [b W];
            upper_a{len+1} = [b W];
            
            % new lower and uppper bounds
            lb = obj.lb;
            ub = obj.ub;
            lb{len+1} = obj.lb_backSub(lower_a, upper_a);
            ub{len+1} = obj.ub_backSub(lower_a, upper_a);
            
            R = RStar0(lower_a, upper_a, lb, ub, obj.iter); 
        end
        
        % lower bound back-substitution
        function lb = lb_backSub(obj, lower_a, upper_a)
            maxIter = obj.iter;
            len = length(upper_a);
            [nL, mL] = size(upper_a{len});
            alpha = upper_a{len}(:,2:end);
            lower_v = zeros(nL, 1);
            upper_v = upper_a{len}(:,1);
            
            % b[s+1] = v' + sum( max(0,w[j]')*lower_a[j] + min(w[j]',0)*upper_a[j}] ) for j is element of k and for k < i
            % iteration until lb' = b[s'] = v''
            len = len - 1;
            iter = 0;
            while (len > 1 && iter < maxIter)
                [nL, mL] = size(upper_a{len});
                dim = nL;
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);

                lower_v = max_a * lower_a{len}(:,1) + lower_v;
                upper_v = min_a * upper_a{len}(:,1) + upper_v;
                
                alpha = max_a * lower_a{len}(:,2:end) + ...
                        min_a * upper_a{len}(:,2:end);

                len = len - 1;
                iter = iter + 1;
            end
            
            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            [lb1,ub1] = getRanges_L(obj,len);
            lb = max_a * lb1 + lower_v + ...
                 min_a * ub1 + upper_v;
        end
        
        % upper bound back-substituion
        function ub = ub_backSub(obj, lower_a, upper_a)
            maxIter = obj.iter;
            len = length(upper_a);
            [nL, mL] = size(upper_a{len});
            alpha = upper_a{len}(:,2:end);
            lower_v = zeros(nL, 1);
            upper_v = upper_a{len}(:,1);
            
            % c[t+1] = v' + sum( max(0,w[j]')*upper_a[j] + min(w[j]',0)*lower_a[j}] ) for j is element of k and for k < i
            % iteration until ub' = c[t'] = v''
            len = len - 1;
            iter = 1;
            while (len > 1 && iter < maxIter)
                [nL, mL] = size(upper_a{len});
                dim = nL;
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);
                
                lower_v = min_a * lower_a{len}(:,1) + lower_v;
                upper_v = max_a * upper_a{len}(:,1) + upper_v;
                
                alpha = min_a * lower_a{len}(:,2:end) + ...
                        max_a * upper_a{len}(:,2:end);

                len = len - 1;
                iter = iter + 1;
            end

            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            [lb1,ub1] = getRanges_L(obj,len);
            ub = min_a * lb1 + lower_v + ...
                 max_a * ub1 + upper_v;
        end
        
        % get the lower and upper bound of a current layer at specific
        % position
        function [lb,ub] = getRange(obj, i)
            if i > obj.dim
                error('i should not exceed dimnesion');
            end
            
            len = length(obj.lb);
            lb = obj.lb{len}(i);
            ub = obj.ub{len}(i);
        end
        
        % get lower and upper bounds of a current layer
        function [lb,ub] = getRanges(obj)
            len = length(obj.lb);
            lb = obj.lb{len};
            ub = obj.ub{len};
        end
        
        % get lower and upper bound of a specific layer
        function [lb,ub] = getRanges_L(obj, len)
            numL = length(obj.lower_a);
            if len > numL
                error('range request should be layers within iteration');
            end
            lb = obj.lb{len};
            ub = obj.ub{len};
        end
        
        % convert to Polyhedron
        function P = toPolyhedron(obj)
            [lb, ub] = getRanges(obj);
            P = Polyhedron('lb', [lb], 'ub', [ub]);
        end
        
        % convert to Zonotope
        function Z = toZono(obj)
            [lb,ub] = getRanges(obj);
            c = (ub + lb)*0.5;
            v = (ub - lb)*0.5;
            dim = obj.dim;
            
            V = [];
            for i=1:dim
                V1 = zeros(dim, 1);
                V1(i) = v(i);
                V = [V V1];   
            end
            
            Z = Zono(c, V);
        end
        
        % convert to Star
        function S = toStar(obj)
            [lb,ub] = getRanges(obj);
            c = (ub + lb)*0.5;
            v = (ub - lb)*0.5;
            dim = obj.dim;
            
            V = [];
            for i=1:dim
                V1 = zeros(dim, 1);
                V1(i) = v(i);
                V = [V V1];   
            end
            C = [eye(dim); -eye(dim)];
            d = ones(2*dim, 1); 
            
            S = Star([c V], C, d);
        end
        
        % plot RStar0 set
        function plot (varargin)
            color = 'red';
            switch nargin
                case 1
                    obj = varargin{1};
                case 2
                    obj = varargin{1};
                    color = varargin{2};
                otherwise
                    error('Invalid number of input arguments, should be 1 or 2');
            end
            
            if strcmp(color, 'rand')
                c_rand = rand(1,3);
            else
                c_rand = color;
            end
            
            P = obj.toPolyhedron;
            plot(P, 'color', c_rand);
        end
        
        function bool = isEmptySet(obj)
            S = obj.toStar;
            bool = S.isEmptySet;
        end
    end
end