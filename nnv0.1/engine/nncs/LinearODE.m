classdef LinearODE
    %Linear system represented by ODE x' = Ax + Bu, y = Cx + Du
    %   Linear ODE system class
    
    properties
        A = []; % system matrix A
        B = []; % system matrix B
        C = []; % system matrix C
        D = []; % system matrix D
        dim = 0; % system dimension
        nI = 0; % number of inputs
        nO = 0; % number of outputs
    end
    
    methods
        
        % constructor
        function obj = LinearODE(A, B, C, D)       
            
            if isempty(A)
                error('Matrix A is empty');
            else 
                [nA, mA] = size(A);
                if nA ~= mA
                    error('A should be a square matrix');
                end
            end
            
            if ~isempty(B)
                [nB, mB] = size(B);
                if nA ~= nB
                    error('A and B are inconsistent');
                end
            end
            
            if ~isempty(C)
                [nC, mC] = size(C);
                if mC ~= nA
                    error('A and C are inconsistent');
                end
            end
            
            if ~isempty(C) && ~isempty(D)
               [nD, ~] = size(D); 

                if nD ~= nC
                    error('C and D are inconsistent');
                end
            end

            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            obj.nI = mB;
            obj.nO = nC;
            obj.dim = nA;
                
        end
        
        % simulate the time response of linearODE system
        function [y, t, x] = simulate(obj, u, t, x0)
            % @x0: initial vector
            % @u : control input
            % @y : output value
            
            
            sys = ss(obj.A, obj.B, obj.C, obj.D);
            

            if size(x0, 1) ~= obj.dim || size(x0, 2) ~= 1
                error('initial vector has an inappropriate dimension');
            end
            
            [y, t, x] = lsim(sys, u, t, x0);
            lsim(sys, u, t, x0); % just for plot      
            
        end
        
        % step response of linear ODE system      
        function [y, t, x] = step(obj, Tfinal)
            % @Tfinal: final time of simulation
            % @y: output 
            % @t: time sequence of simulation
            % @x: states 
            
            sys = ss(obj.A, obj.B, obj.C, obj.D);
            [y, t, x] = step(sys, Tfinal);
            step(sys, Tfinal); % just for plot 
            
        end
        
        % response to intial condition
        function [y, t, x] = initial(obj, x0, Tfinal)
            % @x0: initial condition
            % @Tfinal: final simulation time
            % @y: output
            % @t: time sequence
            % @x: states 
            
            sys = ss(obj.A, [], obj.C, []);
            [y, t, x] = initial(sys, x0, Tfinal);
            initial(sys, x0, Tfinal); % just for plot
            
        end
        
        % convert to discrete ODE       
        function sysd = c2d(obj, Ts)
            % @Ts: sampling time
            % @sysd: discrete linear ODE
            
            sys = ss(obj.A, obj.B, obj.C, obj.D);
            
            sys1 = c2d(sys, Ts); % convert to discrete model
            
            sysd = DLinearODE(sys1.A, sys1.B, sys1.C, sys1.D, Ts);
            
        end
                
    end
    
    methods(Static)
    
        % simulation of dot{x} = Ax using direct method
        function X = simDirect(A, x0, h, N)
            % @A: system matrix
            % @x0: initial states
            % @h: timeStep for simulation
            % @N: number of steps
            % @X : simulation result X = [x0, x1, ..., xN]
            
            % author: Dung Tran 
            % date : May/30/2017
            
            n = size(A,1);

            E = expm(h*A);
            X = zeros(n,N);
            x01 = x0;
            for i=1:N
                X(:,i) = E*x01;
                x01 = X(:,i);
            end
            
            X = [x0 X];
            
        end
        
        % simulation of dot{x} = Ax using Krylov subspace method
        function X = simKrylov(A, x0, h, N, m)
            % @A : the input matrix
            % @x0 : the initial state vector 
            % @h: time Step for simulation
            % @N: number of Step we want to simulate
            % @m : the number of basic vectors for Krylov supspace Km, Vm = [v1, .., vm]
            % @X : simulation result X = [x0, x1, ..., xN]
            
            % author: Dung Tran 
            % date : May/30/2017
            
            function [Q,H] = arnoldi(A,q1,m)
                %ARNOLDI    Arnoldi iteration
                %   [Q,H] = ARNOLDI(A,q1,M) carries out M iterations of the
                %   Arnoldi iteration with N-by-N matrix A and starting vector q1
                %   (which need not have unit 2-norm).  For M < N it produces
                %   an N-by-(M+1) matrix Q with orthonormal columns and an
                %   (M+1)-by-M upper Hessenberg matrix H such that
                %   A*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*E_M',
                %   where E_M is the M'th column of the M-by-M identity matrix.

                n = length(A);
                if nargin < 3, m = n; end
                q1 = q1/norm(q1);
                Q = zeros(n,m); Q(:,1) = q1;
                H = zeros(min(m+1,m),n);

                for k=1:m
                    z = A*Q(:,k);
                    for i=1:k
                        H(i,k) = Q(:,i)'*z;
                        z = z - H(i,k)*Q(:,i);
                    end
                    if k < n
                       H(k+1,k) = norm(z);
                       if H(k+1,k) == 0, return, end
                       Q(:,k+1) = z/H(k+1,k);
                   end
                end
            end


            [mA] = size(A,1);
            % Krylov supspace method
            Im = eye(m);
            e1 = Im(:,1); % 1fst unit vector of R^m 
            X = zeros(mA,N);
            [V,H] = arnoldi(A,x0,m);
            beta  = norm(x0);
            Vm = V(:,1:m);
            Hm = H(1:m,1:m);
            Vmm = beta*Vm;
            Hms = h*Hm;

            for i=1:N
                P = expm(i*Hms);
                X(:,i) = Vmm*P*e1;
            end
            
            X = [x0 X];
                   
        end
        
        % simulation of dot{x} = Ax using ode45
        function X = simOde45(A, x0, h, N)
            % @A: system matrix
            % @x0: initial states
            % @h: timeStep for simulation
            % @N: number of steps
            % @X : simulation result X = [x0, x1, ..., xN]; 
            
            % author: Dung Tran 
            % date : May/30/2017
            
            function dydt = odefun(y,A)
                %This function generate ode function from matrix A
                [mA,nA] = size(A);
                dydt = zeros(mA,1);
                for i=1:mA
                    for j=1:nA
                        dydt(i)= dydt(i) + A(i,j)*y(j);
                    end
                end
            end
            
            timeSpan = 0:h:h*N;
            [~,y] = ode45(@(t,y) odefun(y,A), timeSpan, x0);
            X = y';
            
        end
        
        
        % simulation-based reachability analysis for dot{x} = Ax using
        % direct method and star set
        function R = simReachDirect(A, X0, h, N)
            % @A: system matrix
            % @X0: initial set of states (a Star set)
            % @h: timeStep for simulation
            % @N: number of steps
            % @X : reachable set X = [X0, X1, ..., XN] (Star set)
            
            % author: Dung Tran 
            % date : May/30/2017
            
            if ~isa(X0, 'Star')
                error('Initial set is not a Star');
            end
            
            m = size(X0.V, 2); % number of basic vectors
            Z = cell(1, m);
            
            for i=1:m
                Z{1, i} = LinearODE.simDirect(A, X0.V(:, i), h, N);
            end
           
            R = [];
            for i=1:N+1
                V = [];
                for j=1:m
                    V = [V Z{1, j}(:, i)];                  
                end
                R = [R Star(V, X0.C, X0.d)]; % reachable set
            end
                        
        end
        
        
              
    end
    
end

