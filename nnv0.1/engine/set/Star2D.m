classdef Star2D
    % A generalization of Star class to 2-dimensional space
    % Used to represent attacked images
    % We can generalize Star to arbitrary high dimensional space
    % Dung Tran: 12/17/2018
    
    % ====================================================================%
    %                   Definition of Star2D
    % 
    % A 2D star set S is defined by: 
    % S = {x| x = C + a[1]*V[1] + a[2]*V[2] + ... + a[n]*V[n]
    %           = V * b, V = {c V[1] V[2] ... V[n]}, 
    %                    b = [1 a[1] a[2] ... a[n]]^T                                   
    %                    where C*a <= d, constraints on a[i]}
    % where, C, V[i] are 2D matrices with the same dimension, i.e., 
    % C and V[i] \in R^{m x n}
    % C : is called the center matrix and V[i] is called the basic matrix 
    %
    % The notion of Star2D is more general than the original Star set where
    % the C and V[i] are vectors. 
    % 
    % Dimension of Star2D is the dimension of the center matrix C
    % ====================================================================%
    
    
    properties
        V = []; % cell array of basic matrices of Star2D
        C = []; % constraints matrix of Star2D
        d = []; % constraints vector of Star2D
        dim = []; % dimension of Star2D
        nVar = 0; 
    end
    
    methods
        
        % constructor
        function obj = Star2D(varargin)
            % @V: a cell array of center matrix and basic matrices 
            % @C: a constraints matrix
            % @d: a constraints vector
            
            switch nargin
                case 3
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                                     
                    n = size(V);
                    m = size(C);
                    l = size(d);

                    if n(1) ~= 1 
                        error('cell array of basic matrices should have one row');
                    end
                    if n(2) ~= m(2) + 1
                        error('Inconsistency between the number of basic matrices and the number of constraints varibles');
                    end
                    if m(1) ~= l(1)
                        error('Inconsistency between the constraint matrix and constraint vector');
                    end
                    if l(2) ~= 1
                        error('Constraint vector should have one column');
                    end

                    % checking consistency between center matrix and basic matrices
                    o = size(V{1});

                    for i=1:n(2)
                        oo = size(V{i});                
                        if oo(1) ~= o(1) || oo(2) ~= o(2)
                            error('Inconsistency between center matrix and basic matrices');
                        end
                    end

                    obj.V = V; 
                    obj.C = C;
                    obj.d = d;
                    obj.dim = o; 
                    obj.nVar = m(2); 

                otherwise
                    % create empty Star2D
                    obj.V = [];
                    obj.C = [];
                    obj.d = [];
                    obj.nVar = 0;
            end
            
            
        end
        
        
    end
end

