classdef Reduction
    % Reduction class
    % A class contains methods for reducing number of contraints of
    % polyhedrons or doing quick convex-hull operations of two polyhedrons
    
    properties
    end
    
    methods(Static)
        
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
            
            if n1 <= n2
                V = V1;
                P = Polyhedron(V');
                for i=1:n2
                    if ~P.contains(V2(:, i))
                        V = [V V2(:, i)];
                        P = Polyhedron(V');
                    end
                end
                                
            else
                V = V2;
                P = Polyhedron(V');
                for i=1:n1
                    if ~P.contains(V1(:, i))
                        V = [V V1(:, i)];
                        P = Polyhedron(V');
                    end
                end
            end
            fprintf('\nNumber of contraints = %d', size(P.A, 1));
            
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
            for i=1:l
                if ~I(i).isBounded
                    error('I(%d) is not bounded', i);
                end
                I(i).minVRep();
                V = [V I(i).V'];
            end
            V1 = unique(V', 'rows', 'stable');
            V = V1';
            
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

