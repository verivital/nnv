classdef ODEblockLayer
    %ODEBLOCKLAYER class is a layer class for Verification of NeuralODEs
    % Contains contructor, simulation and reachability analysys methods
    % 
    % Syntax:
    %    layer = ODEblockLayer(varargin)
    % 
    % Example:
    %    layer = ODEblockLayer(odeblock,1)
    %      where:
    %          odeblock = LinearODE(Aout,Bout,Cout,D,tfinal,numSteps);
    %
    % Example2:
    %    layer = ODEblockLayer(odeblock,1, 0.1, false)
    % 
    % Diego Manzanas: February 15th, 2021
    
    properties
        tfinal = 1; % Final integration time (default = 1)
        tstep = 1/10; % Initialize time step used to return time series
        odemodel = []; % ODE reachability model (linear/nonlinear)
        time_series = []; % Boolean: ture or false. Represents a time-series neural ODE or just another layer (output at tf only)
    end
    
    methods % main methods 

        % Constructor of the class
        function obj = ODEblockLayer(varargin)
        % layer = ODEblockLayer(odemodel,tfinal,tstep)
            if nargin == 1
                obj.odemodel = varargin{1};
                obj.tfinal = 1;
                obj.tstep = 1/10;
                obj.time_series = logical(false);
            elseif nargin == 2
                obj.odemodel = varargin{1};
                obj.tfinal = varargin{2};
                obj.tstep = 1/10;
                obj.time_series = logical(false);
            elseif nargin == 3
                obj.odemodel = varargin{1};
                obj.tfinal = varargin{2};
                obj.tstep = varargin{3};
                obj.time_series = logical(false);
            elseif nargin == 4
                obj.odemodel = varargin{1};
                obj.tfinal = varargin{2};
                obj.tstep = varargin{3};
                obj.time_series = logical(varargin{4});
            else
                error('Wrong number of inputs.');
            end
            
            % Assert inputs are correct
            % 1) ODEmodel
            if ~ any(strcmp(class(obj.odemodel),{'LinearODE', 'NonLinearODE', 'LinearODE_cora'}))
                error('The ODE layer must be a continuous-time model.')
            end
            % 2) Check tfinal and tstep data type and correctness
            if ~ all(strcmp('double',{class(obj.tstep),class(obj.tfinal)}))
                error('Time step and final time must be type double')
            elseif obj.tstep > obj.tfinal
                error('Time step must be smaller than time final.')
            elseif rem(obj.tfinal,obj.tstep)
                error('Wrong time step. Final time must be divisible by time step.');
            end
            % 3) Check time_series is of type logical
            if ~ islogical(obj.time_series)
                error('Wrong input type time_series. It must be false or true.')
            end
        end
        
        % Simulate the ODEblock
        function y = evaluate(obj,input)
            x0 = input;
            u = zeros(obj.odemodel.dim,1);
            if ~ isa(obj.odemodel,'LinearODE')
                [~,x] = obj.odemodel.evaluate(x0,u);
                if obj.time_series
                    y = x';
                else
                    y = x(end,:)';
                end
            else
                reachStep = obj.odemodel.controlPeriod/obj.odemodel.numReachSteps;
                tvec = 0:reachStep:obj.odemodel.controlPeriod;
                u = ones(length(tvec),obj.odemodel.nI);
                [x,~] = obj.odemodel.simulate(u,tvec,x0);
                if obj.time_series
                    y = x';
                else
                    y = x(end,:)';
                end
            end
        end
        
        % Compute reachability analysis 
        function Rf = reach(obj,X,varargin)
            if contains(class(X),'Image')
                X = X.toStar;
            elseif ~ contains(class(X),'Star')
                error('Wrong input set representation. Only Star and ImageStar are allowed as inputs')
            end

            % Different functions for different systems
            if isa(obj.odemodel,'LinearODE')
                U = Star(1,1);
                reachStep = obj.odemodel.controlPeriod/obj.odemodel.numReachSteps;
                R = obj.odemodel.simReach('direct',X,U,reachStep,obj.odemodel.numReachSteps);
                if obj.time_series
                    Rf = R;
                else
                    Rf = R(end);
                end
            else
                U = Star(0,0);
                R = obj.odemodel.stepReachStar(X,U);
                if obj.time_series
                    obj.odemodel.get_interval_sets;
                    Rf = obj.odemodel.intermediate_reachSet;
                else
                    Rf = R;
                end
            end
            
        end
    
    end

end

