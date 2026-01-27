classdef Partition
    % This class partition an input box into n^2 boxes
    % A box is defined by B = [lb ub] \in R^{n x 2}, lb, ub is the
    % lowerbound and upperbound vector
    
    % Dung Tran: 9/20/2018
    
    properties
    end
    
    methods(Static)
        
        % step devide will devide one box into two boxes
        function B = stepDivide(I, idx)
            % @I: input boxes
            % idx: the index of the element that we want to divide
            
            if iscell(I)
                n = length(I);
                B = cell(2*n, 1);
                for i=1:n
                    Bi = I{i, 1};
                    lb = Bi(:,1);
                    ub = Bi(:,2);
                    
                    lb1 = lb;
                    ub1 = ub;
                    ub1(idx,1) = 0.5*(lb(idx,1) + ub(idx,1));

                    lb2 = lb;
                    ub2 = ub;
                    lb2(idx,1) = ub1(idx,1);

                    B1 = [lb1 ub1];
                    B2 = [lb2 ub2];
                    B{2*i-1, 1} = B1;
                    B{2*i, 1} = B2;
                end
    
                
            else
                error('Input box should be in a cell');
            end
                        
        end
        
        % partition a box I into a set of smaller boxes
        function P = partition_box(I, N)
            % @I: input polyhedra (ideally a box) if not we get an
            % over-approximation of the polyhedra
            % @N: number of step devide for each element
                        
            
            % @P: a cell of partitioned boxes
            %     |P| = n^(2^N), n is the dimension of the box
            
            
            if ~isa(I, 'Box')
                I.outerApprox;
                n = size(I.A, 2);
                lb = I.Internal.lb;
                ub = I.Internal.ub;
            else
                lb = I.lb;
                ub = I.ub;
                n = size(lb, 1);
            end
            
            B = cell(1,1);
            B{1, 1} = [lb ub];   % put box into a cell
                 
            for i=1:n
                if lb(i) ~= ub(i)
                    for j=1:N
                        B = Partition.stepDivide(B, i);
                    end
                end
                
            end
            
            L = length(B);
            P = [];
           
            for i=1:L
                P = [P Box(B{i,1}(:,1), B{i,1}(:,2))];
            end
            
                
        end
        
        
        % approximate a convex polyhedra by boxes
        function B = approximateByBoxes(I,n, nB, parallel)
            % @I: a polyhedra
            % @n: divide the box bounding I by n times.
            % @N : number of boxes want to use to approximate I
            % @B: an array of boxes
            
            % only work for small-dimensional polyhedron <= 4
            
            B1 = Partition.partition_box(I, n);
            L = length(B1);
            B = [];
            if strcmp(parallel, 'parallel')
                parfor i=1:L
                    if ~isEmptySet(B1(i) & I)
                        B = [B B1(i)];
                    end 
                end
                
                B = Reduction.merge_box(B, nB, parallel);
                
            elseif strcmp(parallel, 'single')
                for i=1:L
                    if ~isEmptySet(B1(i) & I)
                        B = [B B1(i)];
                    end 
                end
                
                B = Reduction.merge_box(B, N, parallel);
                
            else
                error('Unknown parallel computing option');
            end
            
        end
        
        
    end
    
end

