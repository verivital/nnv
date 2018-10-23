classdef DLinearODE
    %Discrete Linear ODE class
    %   represents discrete linear ode system: x[k + 1] = Ax[k] + Bu[k],
    %                                          y[k + 1] = Cx[k] + Du[k]
    %   Dung Tran: 10/21/2018
    
    properties
        A = []; % system matrix A
        B = []; % system matrix B
        C = []; % system matrix C
        D = []; % system matrix D
        nI = 0; % number of inputs
        nO = 0; % number of outputs
        dim = 0; % system dimensions
        Ts = 0 ; % sampling time
    end
    
    methods
        
        % constructor
        function obj = DLinearODE(A, B, C, D, Ts)
            
                
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
            obj.Ts = Ts;
                
            
        end
        
        % simulate the time response of discrete linearODE system
        function [y, t, x] = simulate(obj, u, t, x0)
            % @x0: initial vector
            % @u : control input
            % @y : output value
            
            
            sys = ss(obj.A, obj.B, obj.C, obj.D, obj.Ts);
            

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
            
            sys = ss(obj.A, obj.B, obj.C, obj.D, obj.Ts);
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
            
            sys = ss(obj.A, [], obj.C, [], obj.Ts);
            [y, t, x] = initial(sys, x0, Tfinal);
            initial(sys, x0, Tfinal); % just for plot
            
        end
        
        % convert to continuous linear ODE
        function sysc = d2c(obj)
            
            sysd = ss(obj.A, obj.B, obj.C, obj.D, obj.Ts);
            sys1 = d2c(sysd);
            
            sysc = LinearODE(sys1.A, sys1.B, sys1.C, sys1.D);
            
        end
        
        % reachability analysis of DlinearODE using star set
        function R = stepReachStar(obj, I, U)
            % @I: set of initial condition 
            % @U: set of control input
            % @R: state reachable set (a star set)
            
            if ~isa(I, 'Star')
                error('Set of initial condition is not a star set');
            end
            if ~isa(U, 'Star')
                error('Set of control input is not a star set');
            end
            
            % R = AI + BU
            R1 = I.affineMap(obj.A, []);
            R2 = U.affineMap(obj.B, []);
            
            R = R1.MinkowskiSum(R2);
            
        end
        
        % reachability analysis of DLinearODE using zonotope
        function Z = stepReachZono(obj, I, U)
            % @I: set of initial condition (a zonotope)
            % @U: set of control input (a zonotope)
            % @R: state reachable set (a zonotope)
            
            if ~isa(I, 'Zono')
                error('Set of initial condition is not a zonotope');
            end
            if ~isa(U, 'Zono')
                error('Set of control input is not a zonotope');
            end
            
            Z1 = I.affineMap(obj.A, []);
            Z2 = U.affineMap(obj.B, []);
  
            Z = Z1.MinkowskiSum(Z2);          
        end
        
        
    end
    
end

