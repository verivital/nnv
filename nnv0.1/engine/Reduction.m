classdef Reduction
    % Reduction class
    % A class contains methods for reducing number of contraints of
    % polyhedrons or doing quick convex-hull operations of two polyhedrons
    
    % Dung Tran: 6/9/2018
    
    properties
    end
    
    methods(Static)
        
        % Recursively merge a number of polyhedrons into one with more
        % constraints
        % This function follows the "divide and conquer idea"
        % The complexity of this algorithm is O(nlog(n)) where n is the
        % number of polyhedrons.
        function P = recursiveMerge(I, parallel)
            % @I: an array of input polyhedrons
            % @paprallel: = 'parallel' use parallel computing 
            %             = 'single' use single core for computing
            
            %if length(I) == 1
            %    return I;
            %elseif length(I) == 2
            %    return Reduction.merge(I, parallel);
            %else 
            %    R = Reduction.stepMerge(I, parallel);
            %    return Reduction.recursiveMerge(R, parallel);
            %end
            
            startTime = tic;
            P = I;
            i = 0;
            while length(P) > 1
                i = i + 1;
                fprintf('\nPerforming Step Merging %d (merging %d polyhedrons into %d polyhedrons)...', i, length(P), ceil(length(P)/2));
                P = Reduction.stepMerge(P, parallel);
                
            end
            runtime = toc(startTime);
            fprintf('\nTotal number of step merging = %d', i);
            fprintf('\nTotal merging time = %.5f seconds', runtime);
            fprintf('\nThe merged polyhedron has %d inequality constraints and %d equality constraints', size(P.A, 1), size(P.A, 2));
            
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
            
            % using parallel computing
            if strcmp(parallel, 'parallel')

                if ~isempty(n1_ind)
                    for i=1:length(n1_ind)
                        A = vertcat(A, P1.A(n1_ind(i), :));
                        b = vertcat(b, P1.b(n1_ind(i), :));
                    end
                end


                if ~isempty(m1_ind)
                    for i=1:length(m1_ind)
                        Ae = vertcat(Ae, P1.Ae(m1_ind(i), :));
                        be = vertcat(be, P1.be(m1_ind(i), :));
                    end
                end


                if ~isempty(n2_ind)
                    for i=1:length(n2_ind)
                        A = vertcat(A, P2.A(n2_ind(i), :));
                        b = vertcat(b, P2.b(n2_ind(i), :));
                    end
                end


                if ~isempty(m2_ind)
                    for i=1:length(m2_ind)
                        Ae = vertcat(Ae, P2.Ae(m2_ind(i), :));
                        be = vertcat(be, P2.be(m2_ind(i), :));
                    end
                end

            % using single core    
            elseif strcmp(parallel, 'single')
            
                if ~isempty(n1_ind)
                    for i=1:length(n1_ind)
                        A = vertcat(A, P1.A(n1_ind(i), :));
                        b = vertcat(b, P1.b(n1_ind(i), :));
                    end
                end


                if ~isempty(m1_ind)
                    for i=1:length(m1_ind)
                        Ae = vertcat(Ae, P1.Ae(m1_ind(i), :));
                        be = vertcat(be, P1.be(m1_ind(i), :));
                    end
                end


                if ~isempty(n2_ind)
                    for i=1:length(n2_ind)
                        A = vertcat(A, P2.A(n2_ind(i), :));
                        b = vertcat(b, P2.b(n2_ind(i), :));
                    end
                end


                if ~isempty(m2_ind)
                    for i=1:length(m2_ind)
                        Ae = vertcat(Ae, P2.Ae(m2_ind(i), :));
                        be = vertcat(be, P2.be(m2_ind(i), :));
                    end
                end

                
            else 
                error('Unknown parallel option');
            end
            
            
            A = vertcat(A, B.A);
            b = vertcat(b, B.b);
            Ae = vertcat(Ae, B.Ae);
            be = vertcat(be, B.be);
            
            % the final merged polyhedron
            P = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be).minHRep();
            
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
        function P = fastHull(P1, P2)
            % the convexHull operation in mpt toolbox may crash when
            % the dimensions of P1 and P2 increase
            % This fast hull implement a simple algorithm to construct
            % A convex hull of two polyhedrons using their vertices
            % This function exploits the fact that a bounded polyhedron has
            % no rays.
            % @P1: the first polyhedron
            % @P2: the second polyhedron
            % %P : the hull of two polyhedron
            
            V1 = P1.V';
            V2 = P2.V'; 
            
            [m1, n1] = size(V1); % number of vertices of P1
            [m2, n2] = size(V2); % number of vertices of P2
            
            if m1 ~= m2
                error('Inconsistent dimensions between two polyhedrons');
            end
            
            V = [V1 V2];
            P = Polyhedron(V').minHRep();           
            %fprintf('\nNumber of contraints = %d', size(P.A, 1));
            
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
            
            tic;
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
                    for i=1:n
                        P = I((i-1)*m + 1 : i*m);
                        V = [];
                        for j=1:m
                            V = [V P(j).V'];
                        end
                        R = [R Polyhedron(V').outerApprox];
                    end

                    if n*m < N
                        P = I(n*m:N);
                        V = [];
                        for j=1:N-n*m
                            V = [V P(j).V'];
                        end
                        R = [R Polyhedron(V').outerApprox];
                    end

                else
                    
                    P = I(1:N);
                    V = [];
                    for j=1:N
                        V = [V P(j).V'];
                    end
                    R = [R Polyhedron(V').outerApprox];
                end
                
                    
            end
            t = toc;
                       
            
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
            for i=1:l
                P(i) = I(i);
            end
            U = PolyUnion(P);
            B = U.outerApprox;
                      
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

