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
        function R = stepReach(I, index, xmin, xmax)
            % @I : input set, a polyhedron
            % @i : index of current x[index] of current step
            % @xmin: min of x[index]
            % @xmax: max of x[index]
           
           
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
                b  = vertcat(I.b, [0]);
                R1 = Polyhedron('A', A1, 'b', b, 'Ae', I.Ae, 'be', I.be);
                R2 = Polyhedron('A', A2, 'b', b, 'Ae', I.Ae, 'be', I.be);
                
                Im = eye(dim);
                Im(index, index) = 0;
                R1 = R1.affineMap(Im);
                
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
                        if R1 <= R2
                            R = R2;
                        else
                            R = [R1 R2];
                        end
                    end
                end
                
                
            end
            
        end
        
        % stepReach for multiple Input Sets
        
        function R = stepReachMultipleInputs(I_array, index, xmin, xmax)
            % @I: an array of input sets which are polyhedra
            % @index: index of current x[index] of current step
            % @xmin: min value of x[index]
            % @xmax: max value of x[index]
            % @option: = 'exact' -> compute an exact reach set
            %          = 'approx' -> compute an over-approximate reach set
            
            p = length(I_array);
            R = [];
            
            for i=1:p
                R =[R, ReLU.stepReach(I_array(i), index, xmin, xmax)];                
            end
                     
        end
        
        % stepReach for multiple Input Sets use parallel computing
        
        function R = stepReachMultipleInputs_parallel(I_array, index, xmin, xmax)
            % @I: an array of input sets which are polyhedra
            % @index: index of current x[index] of current step
            % @xmin: min value of x[index]
            % @xmax: max value of x[index]
            % @option: = 'exact' -> compute an exact reach set
            %          = 'approx' -> compute an over-approximate reach set
            
            p = length(I_array);
            R = [];
            
            parfor i=1:p
                R =[R, ReLU.stepReach(I_array(i), index, xmin, xmax)];                
            end            
             
        end
        
        % exact reach of ReLU(x)
        function [R, rn] = reach(I)
            % reachable set computation for ReLU(x)
            % @I: input set which is a polyhedron
            % @R : reachable set of ReLU(x), may be an array of polyhedra
            % @rn: number of ReLU_i (stepReach) operations reduced
            % @option: = 'exact' -> compute an exact reach set
            %          = 'approx' -> compute an over-approximate reach set
            
                        
            I.outerApprox; % find bounds of I state vector
            lb = I.Internal.lb; % min-vec of x vector
            ub = I.Internal.ub; % max-vec of x vector
            map = find(lb < 0); % computation map
            m = size(map, 1); % number of stepReach operations needs to be executed
            rn = size(lb, 1) - m;
            
            In = I;
            for i=1:m
                fprintf('\nPerforming ReLU_%d operation', i);
                In = ReLU.stepReachMultipleInputs(In, map(i), lb(map(i)), ub(map(i)));
            end               
            R = In;
           
        end
        
        % exact reach of ReLU(x) in parallel
        function [R, rn] = reach_parallel(I)
            % reachable set computation for ReLU(x)
            % @I: input set which is a polyhedron
            % @R : reachable set of ReLU(x), may be an array of polyhedra
            % @rn: number of ReLU_i (stepReach) operations reduced
            % @option: = 'exact' -> compute an exact reach set
            %          = 'approx' -> compute an over-approximate reach set
            
                        
            I.outerApprox; % find bounds of I state vector
            lb = I.Internal.lb; % min-vec of x vector
            ub = I.Internal.ub; % max-vec of x vector
            map = find(lb < 0); % computation map
            m = size(map, 1); % number of stepReach operations needs to be executed
            rn = size(lb, 1) - m;
            
            In = I;
            for i=1:m
                fprintf('\nPerforming ReLU_%d operation', i);
                In = ReLU.stepReachMultipleInputs_parallel(In, map(i), lb(map(i)), ub(map(i)));
            end               
            R = In;
           
        end
        
    end
end

