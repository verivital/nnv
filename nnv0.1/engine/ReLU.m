classdef ReLU
    % ReLU operator in NN
    % Dung Tran: 8/21/ 2018
    
    properties
    end
    
    methods(Static)
        
        % evaluation
        function y = evaluate(x)
            n = size(x, 1);
            if size(x, 2) ~= 1
                error('x is not a vector')
            end
            y = zeros(n, 1);
            for i=1:n
                if x(i, 1) < 0
                    y(i, 1) = 0;
                else
                    y(i, 1) = x(i, 1);
                end
            end

        end
        
        % step reach set y = ReLU(x)
        function R = stepReach(I, index, xmin, xmax, option)
            % @I : input set, a polyhedron
            % @i : index of current x[index] of current step
            % @xmin: min of x[index]
            % @xmax: max of x[index]
            % @option: 'exact' -> return an exact reach set
            %           'approx' -> return an over-approximation reach set
            
            dim = size(I.A, 2);  % dimension of input space
            if xmin >= 0
                R = I; 
            elseif xmax < 0 
                Im = eye(dim);
                Im(index, index) = 0;
                R = I.affineMap(Im);    % y[index] = ReLU(x[index]) = 0
            elseif xmin < 0 && xmax >= 0
                Z1 = zeros(1, dim);
                Z1(1, index) = 1;
                Z2 = zeros(1, dim);
                Z2(1, index) = -1;
                
                A1 = vertcat(I.A, Z1);
                A2 = vertcat(I.A, Z2);         
                b = vertcat(I.b, 0);
        
                R1 = Polyhedron('A', A1,'b', b, 'Ae', I.Ae, 'be', I.be); % R1 = I \intersect x[index] <= 0
                R2 = Polyhedron('A', A2, 'b', b, 'Ae', I.Ae, 'be', I.be); % R2 = I \intersect x[index] >= 0
                
                Im = eye(dim);
                Im(index, index) = 0;
                R1 = R1.affineMap(Im);
                
                if R1 <= R2  % if R1 is a subset of R2
                    R = R2;
                elseif R2 <= R1 % if R2 is a subset of R1
                    R = R1;
                else
                    if strcmp(option, 'exact')    
                        R = [R1 R2]; % array of 2 polyhedra
                    elseif strcmp(option, 'approx')
                        if (size(R1.A, 1) >= 15) || (size(R2.A, 1) >= 15)
                            fprintf('Input set dimension is too large for approx option');
                            fprintf('Approx option with compute convex hull of a set of polyhedron for each layer');
                            fprintf('The input set dimension is too large, the convex hull operation may be crashed or very time-consumming');
                            fprintf('Consider choosing approx-box option instead');
                        end
                        P(1) = R1;
                        P(2) = R2;
                        U = PolyUnion('Set', P);
                        R = U.convexHull;
                    else
                        error('Unknown option');
                    end
                end
                  
                
            end
            
        end
        
        % stepReach for multiple Input Sets
        
        function R = stepReachMultipleInputs(I_array, index, xmin, xmax, option)
            % @I: an array of input sets which are polyhedra
            % @index: index of current x[index] of current step
            % @xmin: min value of x[index]
            % @xmax: max value of x[index]
            % @option: = 'exact' -> compute exact reach set
            %          = 'approx' -> compute over-approximate reach set
            
            p = size(I_array, 1);
            R = [];
            for i=1:p
                I = I_array(i);
                O = ReLU.stepReach(I, index, xmin, xmax, option);
                R = [R, O];
            end
 
        end
        
        % reach of ReLU(x)
        function [R, rn] = reach(I, option)
            % reachable set computation for ReLU(x)
            % @I: input set which is a polyhedron
            % @R : reachable set of ReLU(x), may be an array of polyhedra
            % @rn: number of ReLU_i (stepReach) operations reduced
            % @option: = 'exact' -> compute an exact reach set
            %          = 'approx' -> compute an over-approximate reach set
            
            if ~I.isBounded
                error('Input set is not bounded')
            end
            
            if I.isEmptySet
                R = []; % if I is empty, return empty set
                rn = 0;
            else                      
                B = outerApprox(I); % interval hull of the input set
                lb = B.Internal.lb; % min-vec of x vector
                ub = B.Internal.ub; % max-vec of x vector
                map = find(lb < 0); % computation map
                m = size(map, 1); % number of stepReach operations needs to be executed
                rn = size(lb, 1) - m;

                In = I;
                for i=1:m
                    In = ReLU.stepReachMultipleInputs(In, map(i), lb(map(i)), ub(map(i)), option);
                end
                R = In;
            end
        end
        
        
    end
    
end

