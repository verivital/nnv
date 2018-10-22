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
            
            
            if nargin < 3
                error('Missing input matrices to construct a linearODE system');
            elseif nargin == 3
                [nA, mA] = size(A);
                [nB, mB] = size(B);
                [nC, mC] = size(C);


                if nA ~= mA
                    error('A should be a square matrix');
                end
                if nA ~= nB
                    error('A and B are inconsistent');
                end
                if mC ~= nA
                    eror('A and C are inconsistent');
                end

                obj.A = A;
                obj.B = B;
                obj.C = C;
                obj.nI = mB;
                obj.nO = nC;
                obj.dim = nA;

            elseif nargin == 4
                
                [nA, mA] = size(A);
                [nB, mB] = size(B);
                [nC, mC] = size(C);
                [nD, ~] = size(D);
                
                if nA ~= mA
                    error('A should be a square matrix');
                end
                if nA ~= nB
                    error('A and B are inconsistent');
                end
                if mC ~= nA
                    error('A and C are inconsistent');
                end
                                
                if nD ~= nC
                    error('C and D are inconsistent');
                end

                obj.A = A;
                obj.B = B;
                obj.C = C;
                obj.D = D;
                obj.nI = mB;
                obj.nO = nC;
                obj.dim = nA;
                
            elseif nargin > 4
                error('Too many inputs');
            end
                
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
    
end

