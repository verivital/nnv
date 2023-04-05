classdef Reduction
    % Reduction class
    % A class contains methods for reducing number of contraints of
    % polyhedrons or doing quick convex-hull operations of two polyhedrons
    
    % Dung Tran: 6/9/2018
    
    properties
    end
    
    methods(Static)
        
        
        % get vertices from lower-bound and upper-bound vector
        function V = getVertices(lb, ub)
            if size(lb,1) ~= size(ub,1)
                error('Inconsistency between lb and ub');
            end
            
  
            [ms,~] = size(lb);


            % calculate number of vertices in the initial set X0
            pos = zeros(ms,1);
            n = 0;
            for i = 1:ms
                if (lb(i)==ub(i))
                   pos(i) = 0;
                else 
                   pos(i) = i; 
                   n = n+1;        
                end
            end          

            N = 2^n; % number of vertices in the inital set
            T1 = zeros(n,1); 
            k = 0; 
            for i=1:ms
                if(pos(i)~=0)
                    k = k+1; 
                    T1(k) = pos(i);
                end
            end

            T2 = zeros(N,n); 
            for i=0:N-1
                T2(i+1,:) = de2bi(i,n); 
            end

            T3 = zeros(N,n);

            for i=1:N
                for j = 1:n
                    if T2(i,j)== 0
                       T3(i,j) = lb(T1(j));
                    else
                       T3(i,j) = ub(T1(j));
                    end
                end
            end

            X0_set = zeros(ms,N);
            for i = 1:N
                for j = 1:ms
                    if pos(j) == 0
                       X0_set(j,i)= lb(j);
                    else
                       X0_set(j,i) = 0; 
                    end
                end

                for k = 1:n
                    X0_set(T1(k),i) = T3(i,k);
                end

            end

            V = X0_set;             
            
        end
        
        % affineMap this affineMap function in MPT toolbox is not reliable
        % This function use vertices for affineMapping
        
        function P = affineMap(I, W)
            % @I: input polyhedron
            % @W: affine mapping matrix
            
            V = W * I.V';           
            P = Polyhedron('V', V');
        end
        
        
        % find nP smallest polyhedra close to the origin the most
        function P = findSmallestPolyhedraCloseOrigin(I, nP)
            % @I: an array of input polyhedra
            % @nP: the number of polyhedra that is close to the origin
            % @P : the output polyhedra P = {P1, P2} = I, |P2| = nP
            %    : this algorithm divide the input into two sets, the
            %    second set contains nP smalllest polyhedra that is close the most to
            %    the origin
            
            
            % criteria for  clustering: J = 1/2 * d + 1/2 * V 
            
            % d = || (ub + lb)/2||: distance from the center of the box bounding a
            % polyhedron to the origin
            % V = ||ub - lb||, distance between lowerbound and upperbound of the box 
            
            
            n = length(I);
            C = zeros(n, 1);
            for i=1:n
                I(i).outerApprox;
                lb = I(i).Internal.lb;
                ub = I(i).Internal.ub;
                %C(i,1) = norm(ub);
                C(i, 1) = 0.5*norm((ub +lb)/2) + 0.5*norm((ub - lb));
            end

            [~, Index] = sort(C);
            Idx = Index(1:nP); % index of nP smallest number

            P2 = [];
            for i=1:nP
                P2 = [P2 I(Idx(i))];
            end
            
            Idx = Index(nP + 1:n); % index of n-NP largest number
            P1 = [];
            for i=1:length(Idx)
                P1 = [P1 I(Idx(i))];
            end
            
            P = cell(2,1);
            P{1,1} = P1;
            P{2,1} = P2;
            
        end
        
        %clustering a set of polyhedra by dependent variables        
        function P = clusterByDependentVariables(I)
            % @I: an array of input polyhedra
            
            n = length(I);
            m = size(I(1).A, 2); % number of variables (dimension)
            R = zeros(n, m); % matrix to show dependent variables of polyhedra

            for i=1:n

                for j=1:m
                    if sum(I(i).A(:, j) > 0)
                        R(i, j) = 1;
                    end                   
                end

            end

            [C,~,~]=unique(R,'rows','stable');

            P = cell(size(C, 1), 1);

            for i=1:size(C, 1)
                for j=1:n
                    if sum(abs(R(j, :) - C(i, :))) == 0
                        P{i,1} = [P{i, 1} I(j)];
                    end        
                end
            end
                  
        end
        
        
        % clustering a set of polyhedra by overlapness
        function P = clusterByOverlapness(I, nP)
            % @I: an array of input polyhedra
            % @nP: number of groups of the output polyhedra
            
            
            n = length(I);
            B = [];
                       

            for i=1:n
                B = [B I(i).outerApprox];
            end

            m = size(B(1).A, 2);

            C = zeros(n, 2*m);
            for i=1:n
                C(i, :) = [B(i).Internal.lb' B(i).Internal.ub'];
            end

            idx = kmeans(C, nP); % clustering boxes into nP groups

            P = cell(nP, 1);

            for i=1:nP
                for j=1:n
                    if idx(j) == i
                        P{i, 1} = [P{i, 1} I(j)];
                    end
                end
            end

                                      
        end
        
        
        % limit constraint by angle 
        % reference: PhaVer: Algorithmic verification of hybrid system past
        % Hytech. Goran Frehse, STTT 2008
        
        function P = reduceConstraints(I, nc, parallel)
            % @I: an array of bounded polyhedra
            % @nc: number of Constraints of new polyhedra
            % @P: an array of new bounded polyhedra with smaller number of
            % constraints (ideally: nC constraints)
            
            % @parallel: = 'parallel' using parallel computing
            %            = 'single' use single core for computing
            
           
           L = length(I);
           
           P = [];
           
           if strcmp(parallel, 'parallel')
               
               parfor k=1:L
               
               
                   n = size(I(k).A, 1);
                   
                   if n <= nc
                       R = I(k);
                   else
                       
                       C = zeros(n, n);

                       for i=1:n
                           for j=1:n
                               if j~=i
                                   C(i, j) = abs(I(k).A(i, :)* (I(k).A(j,:)'));
                               end

                           end
                       end

                       [row, col] = find(C==max(max(C, [], 2)));

                       A = vertcat(I(k).A(row(1), :), I(k).A(col(1), :));
                       b = vertcat(I(k).b(row(1), :), I(k).b(col(1), :));

                       R = Polyhedron('A', A, 'b', b, 'Ae', I(k).Ae, 'be', I(k).be);

                       maps = zeros(2, 1);
                       maps(1) = row(1);
                       maps(2) = col(1);

                       C(maps(1), :) = 0;
                       C(maps(2), :) = 0;


                       while (size(R.A, 1) < nc || ~R.isBounded || R.isEmptySet) && length(maps) < n

                           m = length(maps);
                           idx_vec = zeros(m, 1);
                           max_val_vec = zeros(m, 1);
                           for i=1:m
                               [max_val_vec(i), idx_vec(i)] = max(C(:, maps(i)));
                           end

                           [max_val, j] = max(max_val_vec);
                           if max_val > 0
                               idx = idx_vec(j);
                           else
                               [row, ~] = find(C==max(max(C, [], 2)));   
                               idx = row(1);
                           end

                           C(idx, :) = 0;
                           maps = vertcat(maps, idx);
                           A = vertcat(R.A, I(k).A(idx, :));
                           b = vertcat(R.b, I(k).b(idx, :));
                           R = Polyhedron('A', A, 'b', b, 'Ae', I(k).Ae, 'be', I(k).be);
                       
                        end

                        
                   end

                    
                    if ~(I(k) <= R)
                        error('reduceConstraints operation gives an error')
                    end

                    P = [P R.minHRep()];

               end
               
           elseif strcmp(parallel, 'single')
               
               for k=1:L
               
               
                    n = size(I(k).A, 1);
                    if n <= nc
                       R = I(k);
                    else
                                               
                       C = zeros(n, n);

                       for i=1:n
                           for j=1:n
                               if j~=i
                                   C(i, j) = abs(I(k).A(i, :)* (I(k).A(j,:)'));
                               end

                           end
                       end

                       [row, col] = find(C==max(max(C, [], 2)));

                       A = vertcat(I(k).A(row(1), :), I(k).A(col(1), :));
                       b = vertcat(I(k).b(row(1), :), I(k).b(col(1), :));

                       R = Polyhedron('A', A, 'b', b, 'Ae', I(k).Ae, 'be', I(k).be);

                       maps = zeros(2, 1);
                       maps(1) = row(1);
                       maps(2) = col(1);

                       C(maps(1), :) = 0;
                       C(maps(2), :) = 0;


                       while (size(R.A, 1) < nc || ~R.isBounded || R.isEmptySet) && length(maps) < n

                           m = length(maps);
                           idx_vec = zeros(m, 1);
                           max_val_vec = zeros(m, 1);
                           for i=1:m
                               [max_val_vec(i), idx_vec(i)] = max(C(:, maps(i)));
                           end

                           [max_val, j] = max(max_val_vec);
                           if max_val > 0
                               idx = idx_vec(j);
                           else
                               [row, ~] = find(C==max(max(C, [], 2)));   
                               idx = row(1);
                           end

                           C(idx, :) = 0;
                           maps = vertcat(maps, idx);
                           A = vertcat(R.A, I(k).A(idx, :));
                           b = vertcat(R.b, I(k).b(idx, :));
                           R = Polyhedron('A', A, 'b', b, 'Ae', I(k).Ae, 'be', I(k).be);
                       
                        end

                        
                    end

                   
                    if ~(I(k) <= R)
                        error('reduceConstraints operation gives an error');
                    end
                    P = [P R.minHRep()];
               end
           else 
               error('Unknown parallel computing option');
           end       
            
        end
        
        % merge polyhedra using boxes        
        function P = merge_box(I, nP, parallel)
            % @I: array of polyhedra
            % @nP: number of polyhera of the output P
            % @parallel: = 'parallel' use parallel computing
            %            = 'single' use single core for computing
            
            n = length(I);
            B = [];
            
            if strcmp(parallel, 'single')
                
                for i=1:n
                    B = [B I(i).outerApprox];
                end

                m = size(B(1).A, 2);

                C = zeros(n, 2*m);
                for i=1:n
                    C(i, :) = [B(i).Internal.lb' B(i).Internal.ub'];
                end

                idx = kmeans(C, nP); % clustering boxes into nP groups

                R = cell(nP, 1);

                for i=1:nP
                    for j=1:n
                        if idx(j) == i
                            R{i, 1} = [R{i, 1} B(j)];
                        end
                    end
                end

                P = [];
                for i=1:nP
                    B = Reduction.hypercubeHull(R{i, 1}); % return a box                    
                    P = [P B];
                end

            elseif strcmp(parallel, 'parallel')
                
                parfor i=1:n
                    B = [B I(i).outerApprox];
                end

                m = size(B(1).A, 2);
                C = zeros(n, 2*m);
                
                for i=1:n
                    C(i, :) = [B(i).Internal.lb' B(i).Internal.ub'];
                end
                
                idx = kmeans(C, nP);
                R = cell(nP, 1);

                for i=1:nP
                    for j=1:n
                        if idx(j) == i
                            R{i, 1} = [R{i, 1} B(j)];
                        end
                    end
                end

                P = [];
                parfor i=1:nP                   
                    B = Reduction.hypercubeHull(R{i, 1}); % return a box                    
                    P = [P B];
                end

            end

            
                  
        end
         
        % merge boxes using box
        
        function P = merge_box2(I, nP, parallel)
            % @I: array of boxes
            % @nP: number of polyhera of the output P
            % @parallel: = 'parallel' use parallel computing
            %            = 'single' use single core for computing
            
            n = length(I);
            B = [];
            
            if strcmp(parallel, 'single')
                
                
                m = length(I(1).lb);
                C = zeros(n, 2*m);
                for i=1:n
                    C(i, :) = [I(i).lb' I(i).ub'];
                end

                idx = kmeans(C, nP); % clustering boxes into nP groups

                R = cell(nP, 1);

                for i=1:nP
                    for j=1:n
                        if idx(j) == i
                            R{i, 1} = [R{i, 1} I(j)];
                        end
                    end
                end

                P = [];
                for i=1:nP
                    P = [P Box.boxHull(R{i, 1})];
                end

            elseif strcmp(parallel, 'parallel')


                m = size(I(1).lb, 1);
                C = zeros(n, 2*m);
                
                for i=1:n
                    C(i, :) = [I(i).lb' I(i).ub'];
                end
                
                idx = kmeans(C, nP);
                R = cell(nP, 1);

                for i=1:nP
                    for j=1:n
                        if idx(j) == i
                            R{i, 1} = [R{i, 1} I(j)];
                        end
                    end
                end

                P = [];
                parfor i=1:nP
                    P = [P Box.boxHull(R{i, 1})];
                end

            end
            
        end
        
        
        % delete redundant polyhedra
        function P = minPolyhedra(I)
            % delete redundant polyhedra or polyhedra that are subsets of
            % others
            % @I: array of polyhedra
            
            L = length(I);
            C = zeros(L, L);

            for i=1:L

                for j=1:L

                    if (I(i) <= I(j))
                        C(i, j) = 1;
                    end

                end
            end

            P = [];

            for i=1:L
                if sum(C(i, :)) == 1
                    P = [P I(i)];
                end
            end
                        
        end
        
        % Recursively merge a number of polyhedrons into one or several polyhedrons with more
        % constraints
        % This function follows the "divide and conquer idea"
        
        function P = recursiveMerge(I, nP, parallel)
            % @I: an array of input polyhedrons
            % @nP: number of polyhedrons of the output 
            %    : nP is usually chosen as the number of cores using in
            %    parallel computing
            % @paprallel: = 'parallel' use parallel computing 
            %             = 'single' use single core for computing
            
           
            startTime = tic; 
            P = I;
            i = 0;
            while length(P) > 2*nP - 1
                i = i + 1;
                fprintf('\nPerforming Step Merging %d (merging %d polyhedrons into %d polyhedrons)...', i, length(P), ceil(length(P)/2));
                P = Reduction.stepMerge(P, parallel);
                
            end
            
            N = length(P);
            fprintf('\nMerging %d polyhedrons into %d polyhedrons', N, nP);
            R1 = [];
            R2 = [];
            for i=1:N
                if i <= nP-1
                    R1 = [R1 P(i)];
                else
                   R2 = [R2 P(i)];
                end
            end
            
            while length(R2) > 1
                R2 = Reduction.stepMerge(R2, parallel);
            end
            
            P = [R1 R2];
            
            runtime = toc(startTime);
            fprintf('\nTotal merging time = %.5f seconds', runtime);
            
            
        end
        
        
        % reducing a half the number of input polyhedrons 
        function P = stepMerge(I, parallel)
            % @I: the input polyhedrons
            % @paprallel: = 'parallel' use parallel computing 
            %             = 'single' use single core for computing
            
            n = length(I); % number of polyhedrons
            
            m = floor(n/2);
            
            P = [];
                
            for i=1:m
                P = [P Reduction.merge(I(2*i - 1), I(2*i), parallel)];
            end
            
            if 2*m < n
                P = [P I(n)];
            end
            
                   
        end
        
        
        % merge two polyhedrons into one polyhedron - more efficient
        % than convex-hull approach. Only use H-representation of
        % polyhedron, this approach is computationally cheap and less
        % conservative than using zonotope (AI^2 toolbox).
        
        function P = merge(P1, P2, parallel)
            % @P1: the first input polyhedron
            % @P2: the second input polyhedron
            % @parallel: = 'parallel' -> using parallel computing
            %            = 'single'  -> using single core for computing
            
            %tic;
            P1.outerApprox;
            P2.outerApprox;
            lb1 = P1.Internal.lb;
            ub1 = P1.Internal.ub;
            lb2 = P2.Internal.lb;
            ub2 = P2.Internal.ub;
            
            m1 = length(lb1);
            m2 = length(lb2);
            
            if m1 ~= m2
                error('Dimension mismatch between two polyhedrons');
            end
            
            lb = zeros(m1, 1);
            ub = zeros(m1, 1);
            
            for i=1:m1
                lb(i) = min(lb1(i), lb2(i));
                ub(i) = max(ub1(i), ub2(i));
            end
            
            B = Polyhedron('lb', lb, 'ub', ub); % box contains two polyhedrons
            
            %t = toc;
            
            %fprintf('\nTime for building an initial box B is %.5f seconds', t);
            
            %tic; 
            n1_ind = [];
            n2_ind = [];
            m1_ind = [];
            m2_ind = [];
            
            n1 = size(P1.A, 1); % number of linear inequalities in P1
            n2 = size(P2.A, 1); % number of linear inequalities in P2
            m1 = size(P1.Ae, 1); % number of linear equalities in P1
            m2 = size(P2.Ae, 1); % number of linear equalities in P2
            
            
            % check which constraints can be added to refine the box B
            
            if strcmp(parallel, 'parallel') % using parallel computing
                
                parfor i=1:n1
                    A = vertcat(P1.A(i, :), B.A);
                    b = vertcat(P1.b(i, :), B.b);
                    B1 = Polyhedron('A', A, 'b', b, 'Ae', B.Ae, 'Be', B.be);

                    if P2 <= B1
                        n1_ind = [n1_ind i];
                    end
                
                end
                
            elseif strcmp(parallel, 'single') % using single core
                for i=1:n1
                    A = vertcat(P1.A(i, :), B.A);
                    b = vertcat(P1.b(i, :), B.b);
                    B1 = Polyhedron('A', A, 'b', b, 'Ae', B.Ae, 'Be', B.be);

                    if P2 <= B1
                        n1_ind = [n1_ind i];
                    end                  
                
                end
            else
                error('Unknown parallel option');
            end
            
                      
            if strcmp(parallel, 'parallel') % using parallel computing
                
                parfor i=1:m1
                    Ae = vertcat(P1.Ae(i, :), B.Ae);
                    be = vertcat(P1.be(i, :), B.be);
                    B1 = Polyhedron('A', B.A, 'b', B.b, 'Ae', Ae, 'Be', be);

                    if P2 <= B1
                        m1_ind = [m1_ind i];
                    end
                end
                
            elseif strcmp(parallel, 'single') % using single core
                
                for i=1:m1
                    Ae = vertcat(P1.Ae(i, :), B.Ae);
                    be = vertcat(P1.be(i, :), B.be);
                    B1 = Polyhedron('A', B.A, 'b', B.b, 'Ae', Ae, 'Be', be);

                    if P2 <= B1
                        m1_ind = [m1_ind i];
                    end
                    
                end
                
            else
                error('Unknown parallel option');
            end
                        
                       
            if strcmp(parallel, 'parallel') % using parallel computing
                
                parfor i=1:n2
                    A = vertcat(P2.A(i, :), B.A);
                    b = vertcat(P2.b(i, :), B.b);
                    B1 = Polyhedron('A', A, 'b', b, 'Ae', B.Ae, 'Be', B.be);

                    if P1 <= B1
                        n2_ind = [n2_ind i];
                    end
                
                end
                
            elseif strcmp(parallel, 'single') % using single core
                
                for i=1:n2
                    A = vertcat(P2.A(i, :), B.A);
                    b = vertcat(P2.b(i, :), B.b);
                    B1 = Polyhedron('A', A, 'b', b, 'Ae', B.Ae, 'Be', B.be);

                    if P1 <= B1
                        n2_ind = [n2_ind i];
                    end
                
                end
            else
                error('Unknown parallel option');
            end   
           
            
            if strcmp(parallel, 'parallel') % using parallel computing
                
                parfor i=1:m2
                    Ae = vertcat(P2.Ae(i, :), B.Ae);
                    be = vertcat(P2.be(i, :), B.be);
                    B1 = Polyhedron('A', B.A, 'b', B.b, 'Ae', Ae, 'Be', be);

                    if P1 <= B1
                        m2_ind = [m2_ind i];
                    end
                end
                
            elseif strcmp(parallel, 'single') % using single core
                
                for i=1:m2
                    Ae = vertcat(P2.Ae(i, :), B.Ae);
                    be = vertcat(P2.be(i, :), B.be);
                    B1 = Polyhedron('A', B.A, 'b', B.b, 'Ae', Ae, 'Be', be);

                    if P1 <= B1
                        m2_ind = [m2_ind i];
                    end
                end
            else
                error('Unknown parallel option');
            end
            
            
            A = [];
            b = []; 
            Ae = [];
            be = [];
            
            
            n = length(n1_ind) + length(n2_ind);
            m = length(m1_ind) + length(m2_ind);
            
            if n >= 1
                for i=1:n
                    
                    if i <= length(n1_ind)
                        
                        A = vertcat(A, P1.A(n1_ind(i), :));
                        b = vertcat(b, P1.b(n1_ind(i), :));
                        A = vertcat(A, B.A);
                        b = vertcat(b, B.b);
                        B1 = Polyhedron('A', A, 'b', b, 'Ae', B.Ae, 'be', B.be);
                        if P1 <= B1 && P2 <= B1
                            B = B1;
                        end
                        A = [];
                        b = [];
                        
                    else
                        
                        A = vertcat(A, P2.A(n2_ind(i - length(n1_ind)), :));
                        b = vertcat(b, P2.b(n2_ind(i - length(n1_ind)), :));
                        A = vertcat(A, B.A);
                        b = vertcat(b, B.b);
                        B1 = Polyhedron('A', A, 'b', b, 'Ae', B.Ae, 'be', B.be);
                        if P1 <= B1 && P2 <= B1
                            B = B1;
                        end
                        A = [];
                        b = [];
                        
                    end
                    
                end
                
            end
            
            if m >= 1
                for i=1:m
                    
                    if i <= length(m1_ind)
                        
                        Ae = vertcat(A, P1.Ae(m1_ind(i), :));
                        be = vertcat(b, P1.be(m1_ind(i), :));
                        Ae = vertcat(Ae, B.Ae);
                        be = vertcat(be, B.be);
                        B1 = Polyhedron('A', B.A, 'b', B.b, 'Ae', Ae, 'be', be);
                        if P1 <= B1 && P2 <= B1
                            B = B1;
                        end
                        Ae = [];
                        be = [];
                        
                    else
                        
                        Ae = vertcat(Ae, P2.Ae(m2_ind(i - length(m1_ind)), :));
                        be = vertcat(be, P2.be(m2_ind(i - length(m1_ind)), :));
                        Ae = vertcat(Ae, B.Ae);
                        be = vertcat(be, B.be);
                        B1 = Polyhedron('A', B.A, 'b', B.b, 'Ae', Ae, 'be', be);
                        if P1 <= B1 && P2 <= B1
                            B = B1;
                        end
                        Ae = [];
                        be = [];
                        
                    end
                    
                end
                
            end
            
            
                       
            % the final merged polyhedron
            P = B.minHRep();
                        
            if ~(P1 <= P) || ~(P2 <= P)
                error('merging operation gives an error');
            end
            
            %t = toc; 
            
            %fprintf('\nTime for refining box B with new constraints to produce a polyhedron P = %.5f seconds', t);
            
            %fprintf('\nNumber of inequalities of the first polyhedron P1 is %d', size(P1.A, 1));
            %fprintf('\nNumber of inequalities of the second polyhedron P2 is %d', size(P2.A, 1));
            %fprintf('\nNumber of inequalities of P1 added to box B = %d', length(n1_ind));
            %fprintf('\nNumber of inequalities of P2 added to box B = %d', length(n2_ind));
            %fprintf('\nNumber of inequalities of the polyhedron P is %d', size(P.A, 1));
            
            %fprintf('\nNumber of equalities of P1 = %d', size(P1.Ae, 1));
            %fprintf('\nNumber of equalities of P2 = %d', size(P2.Ae, 1));
            %fprintf('\nNumber of equalities of P1 added to box B = %d', length(m1_ind));
            %fprintf('\nNumber of equalities of P2 added to box B = %d', length(m2_ind));
            %fprintf('\nNumber of equalities of the polyhedron P = %d', size(P.Ae, 1));
            
        end
        
        % fast hull of two polyhedrons
        function P = fastHull(I)
            % the convexHull operation in mpt toolbox may crash when
            % the dimensions of P1 and P2 increase
            % This fast hull implement a simple algorithm to construct
            % A convex hull of two polyhedrons using their vertices
            % This function exploits the fact that a bounded polyhedron has
            % no rays.
            % @I: the array of polyhedra
            % @P: hull of polyhedra
            
            n = length(I);
            for i=1:n
                V(i) = I(i);
            end
            
            U = PolyUnion('Set', V);
            P = U.convexHull;
                        
        end
        
        % reducing combines with batching technique
        function [R, t] = reduce_and_batch(I, n)
            % batching the input polyhedrons array into two arrays
            % R1 contains all polyhedrons having <= 2*dim inequalities 
            % R2 contains all polyhedrons having > 2*dim inequalities
            
            % Then we batch R2 into n boxes 
            % The array is R = [R1 R2]
            % R contains all polyhedrons having <= 2*dim inequalities ->
            % this will reducing the time for computing reachable set
            % The number of polyhedrons reduced is : N[I] - N[R1] - n
            
            N = length(I);
            R1 = [];
            R2 = [];
            
            dim = size(I(1).A, 2);
            
            tic;
            fprintf('\nDecomposing input set...')
            for i=1:N
                if size(I(i).A, 1) <= 2*dim
                    R1 = [R1 I(i)];
                else
                    R2 = [R2 I(i)];
                end
            end
            t = toc;
            fprintf('\nFinish Decomposition in %.5f seconds', t); 
            
            if ~isempty(R2)
                fprintf('\nBatching input set with > %d inequalities using boxes ...', 2*dim);
                [R2, t] = Reduction.batch(R2, n);
                fprintf('\nFinish batching in %.5f seconnds', t);
                R = [R1 R2];                
            else
                R = R1;
            end
            fprintf('\nNumber of polyhedrons reduced is %d', N - length(R));
                       
        end
        
        % batching polyhedron array
        function [R, t] = batch(I, n)
            % batching polyhedron array I with N polyhedrons
            % into array R containing n polyhedrons, n << N
            % @I : array of polyhedrons 
            % @n : number of batches 
            
            t1 = tic;
            if n < 1
                error('\nNumber of batch should be >= 1');
            end
            N = length(I);
            R = [];
            m = floor(N/n); % number of polyhedrons merged in one batch
       
            if n==1
                R = I;            
            else
                if m >=1
                    for i=1:n - 1
                        P = I((i-1)*m + 1 : i*m);
                        V = [];
                        for j=1:m
                            V = [V P(j).V'];
                        end
                        R = [R Polyhedron(V')];
                    end

                    
                    P = I((n-1)*m + 1:N);
                    V = [];
                    for j=1:N-(n-1)*m
                        V = [V P(j).V'];
                    end
                    R = [R Polyhedron(V')];


                else
                    
                    P = I(1:N);
                    V = [];
                    for j=1:N
                        V = [V P(j).V'];
                    end
                    R = [R Polyhedron(V')];
                end
                
                    
            end
            t = toc(t1);
                       
            
        end
        
        
        % Oriented Rectangular Hulls
        function P = orientedRectangularHull(I)
            % a fast convex-hull algorithm using oriented rectangular
            % Reference: 1) Efficient Representation and Computation of
            % Reachable Sets for Hybrid Systems, Olaf Stursberg, HSCC 2003
            % @I: an array of polyhedrons
            % @P : the oriented rectangular hull of V
            %      @P has 2*n linear constraints, n is dimension of state
            %      vector.
            
            l = length(I);
            V = [];
            %tic;
            for i=1:l
                if ~I(i).isBounded
                    error('I(%d) is not bounded', i);
                end
                V = [V I(i).V'];
            end
            %t = toc;
            %fprintf('\nTime for finding vertices is %.5f', t);
            
            %tic; 
            [n, m] = size(V);            
            x_m = mean(V, 2);   % mean vector (point) of all points
            Vb = zeros(n, m);
            
            if m == 1
                P = Polyhedron(V');
            else
                
                for i=1:m
                    Vb(:, i) = V(:, i) - x_m;
                end
            
                Cov = (Vb * Vb')/(m - 1); % covariance matrix
                [U, S, R] = svd(Cov); % sigular value decoposition of covariance matrix
                               
                A = vertcat(U', -U'); 
                bmax = zeros(n, 1);
                bmin = zeros(n, 1);
                
                for i=1:n
                    L = U(:, i)' * Vb;

                    bmax(i) = max(L) + U(:, i)' * x_m;
                    bmin(i) = -min(L) - U(:, i)' * x_m;
                    
                end
                
                b = vertcat(bmax, bmin);
                
                P = Polyhedron(A, b);
                                
            end
            %t = toc;
            %fprintf('\nTime for finding oriented rectangular hull is %.5f', t);
            
            for i=1:l
                if ~(I(i) <= P)
                    error('Error in computing oriented rectangular hull, use another option, e.g., approx-box');
                end
            end
                   
        end
        
        % hypercubeHull
        
        function B = hypercubeHull(I)
            % a fast convex-hull using hypercube
            % @I: an array of polyhedrons
            % @B: a hypercube Hull of I
            
            l = length(I);
            lb = [];
            ub = [];
            
            for i=1:l
                if ~isa(I(i), 'Box')
                    I(i).outerApprox;
                    lb = [lb I(i).Internal.lb];
                    ub = [ub I(i).Internal.ub];
                else
                    lb = [lb I(i).lb];
                    ub = [ub I(i).ub];
                end
            end
            
            lb = min(lb, [], 2);
            ub = max(ub, [], 2);
            
             B = Box(lb, ub);
        end
        
        
        % zonotope convex hull
        function Z = zonotopeHull(P, V)
            % construct a zonotope enclosing all vertices V
            % the zonotope has a number of generator = n
            % reference: 1) Zonotope as Bounding Volumes
            % @P: set of extreme points P = {p1, p2, ..., pn}
            % @V: a set of unit vectors V = {v1, v2, ..., vk} 
            % @Z: a zonotope Z = (p, <c1v1, c2v2, ..., ckvk>)
            
            [m1, n] = size(P); 
            [m2, k] = size(V);
            
            if m1 ~= m2
                error('Inconsistent dimension between state vector p and generator v');
            else
                m = m1;
            end
            
            % state vector for linear programming
            % x = [c1, ..., ck, b11, b21, bk1, b12, b22, bk2, ..., b1n,
            % ...bkn, y1, y2, ..., ym]'
            % we need to compute p = [y1; y2; ..; ym] and c = [c1; c2; ...; ck]
            % x = [c; b; p], b = [b11; b21; ...; bk1; b12; ...; bk2; b1n; ...; bkn]
            
            % constraints matrix Ax <= b 
            % -ci <= bij <= ci, 1 <= i<= k, 1 <= j <= n -> we have 2*k*n
            % constraints
            
            
            C1 = vertcat(-eye(k), -eye(k));
            C2 = vertcat(-eye(k), eye(k));
            C = C1;
            for i=2:n
                C = vertcat(C, C1);
            end

            D = C2;
            for i=2:n
                D = blkdiag(D, C2);
            end

            A = horzcat(C, D);
            A = horzcat(A, zeros(2 * k * n, m));
            b = zeros(2 * k * n, 1);
            
            % equation matrix Aeq*x = beq;
            
            V1 = V;
            I = eye(m);
            I1 = I;
            for i=2:n
                V1 = blkdiag(V1, V);
                I1 = vertcat(I1, I);
            end
            Aeq = horzcat(V1, I1);
            Aeq = horzcat(zeros(n * m, k), Aeq);
            beq = P(:);
            
            f = zeros(1, k * (n + 1) + m);
            for i=1:k
                f(i) = 1;
            end
  
            x_opt = linprog(f, A, b, Aeq, beq);
            
            p = x_opt(k * (n + 1) + 1: k * (n + 1) + m, 1); % center vector of the zonotope
            c = x_opt(1: k, 1); % ci vector
            G = zeros(m, k);
            for i=1:k
                G(:, i) = c(i) * V(:, i); % generator matrix of the zonotope
            end
            G1 = unique(G', 'rows', 'stable');
            G = G1'; 
            Z = Zonotope(p, G);
            
        end
            
         
        % zonotope hull with random generator
        
        function Z = zonotopeHull_Random(P, k)
            % construct a zonotope enclosing all vertices V
            % the zonotope has a number of generator = n
            % reference: 1) Zonotope as Bounding Volumes
            % @P: set of extreme points P = {p1, p2, ..., pn}
            % @k: a number of random generator 
            % @Z: a zonotope Z = (p, <c1v1, c2v2, ..., ckvk>)
            
            [m, n] = size(P); 
            
            V = rand(m, k);
            for i=1:k
                V(:, i) = V(:, i) / norm(V(:, i));  % unit generator
            end
            
            % state vector for linear programming
            % x = [c1, ..., ck, b11, b21, bk1, b12, b22, bk2, ..., b1n,
            % ...bkn, y1, y2, ..., ym]'
            % we need to compute p = [y1; y2; ..; ym] and c = [c1; c2; ...; ck]
            % x = [c; b; p], b = [b11; b21; ...; bk1; b12; ...; bk2; b1n; ...; bkn]
            
            % constraints matrix Ax <= b 
            % -ci <= bij <= ci, 1 <= i<= k, 1 <= j <= n -> we have 2*k*n
            % constraints
            
            
            C1 = vertcat(-eye(k), -eye(k));
            C2 = vertcat(-eye(k), eye(k));
            C = C1;
            for i=2:n
                C = vertcat(C, C1);
            end

            D = C2;
            for i=2:n
                D = blkdiag(D, C2);
            end

            A = horzcat(C, D);
            A = horzcat(A, zeros(2 * k * n, m));
            b = zeros(2 * k * n, 1);
            
            % equation matrix Aeq*x = beq;
            
            V1 = V;
            I = eye(m);
            I1 = I;
            for i=2:n
                V1 = blkdiag(V1, V);
                I1 = vertcat(I1, I);
            end
            Aeq = horzcat(V1, I1);
            Aeq = horzcat(zeros(n * m, k), Aeq);
            beq = P(:);
            
            f = zeros(1, k * (n + 1) + m);
            for i=1:k
                f(i) = 1;
            end
  
            %options = optimoptions('linprog','Algorithm','interior-point');
            %[x_opt, fval, exitflag, output] = linprog(f, A, b, Aeq, beq, [], [], options);
            [x_opt, fval, exitflag, output] = linprog(f, A, b, Aeq, beq);
            if exitflag == 1
                p = x_opt(k * (n + 1) + 1: k * (n + 1) + m, 1); % center vector of the zonotope
                c = x_opt(1: k, 1); % ci vector
                G = zeros(m, k);
                for i=1:k
                    G(:, i) = c(i) * V(:, i); % generator matrix of the zonotope
                end
                G1 = unique(G', 'rows', 'stable');
                G = G1'; 
                Z = Zonotope(p, G);
            else
                error('The optimization problem is infeasible');
            end
            
            
        end
        
        
    end
    
end

