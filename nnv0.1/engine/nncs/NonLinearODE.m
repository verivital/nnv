classdef NonLinearODE < handle
    % Nonlinear ODE class is a grapper of nonlinearSys class in CORA
    %   Dung Tran: 11/19/2018
    
    properties
        options = []; % option for recahable set computation
        dynamics_func = []; % function to describe dynamics of the system
        dim = 0; % dimension of the system (number of state variable)
        nI = 0; % number of control input
    end
    
    methods
        
        % constructor
        function obj = nonLinearODE(dim, nI, dynamics_func)
            % construct nonLinearODE plant
            % @dim: dimension of the plant
            % @nI: number of input
            % @dynamics_func: dynamics of the plant
            
            % author: Dung Tran
            % date: 11/19/2018
            
            % Note: we construct a nonLinearODE plant with default option
            
            if dim < 1
                error('Dimension should be >=1');
            end
            
            if nI < 0
                error('Number of inputs should be >= 0');
            end
            
            obj.dim = dim;
            obj.nI = nI;
            obj.dynamics_func = dynamics_func;
            
            % default option
            option.tStart = 0;
            option.tFinal = 1; % we provide method to change the option
            
            option.x0 = [];
            option.R0 = []; % initial set for reachability analysis
            option.timeStep = 0.1; % time step for reachable set computation
            option.taylorTerms = 4; % number of taylor terms for reachable sets
            option.zonotopOrder = 2; % zonotope order
            option.intermediateOrder = 5; 
            option.reductionTechnique = 'girard';
            option.errorOrder = 1;
            option.polytopeOrder = 2; % polytope order
            option.reductionInterval = 1e3;
            option.maxError = 0.5*ones(dim, 1);
            option.advancedLinErrorComp = 0;
            option.tensorOrder=2;
            
            obj.options = option; % default option
            
        end
        
        % set taylor terms
        function set_taylorTerms(obj, taylorTerms)
            
            if taylorTerms < 1
                error('Invalid Taylor Terms');
            end
            obj.options.taylorTerms = taylorTerms;
        end
        
        % set zonotope order
        function set_zonotopeOrder(obj, zonotopeOrder)
            if zonotopeOrder < 1
                error('Invalid zonotope order');
            end
            
            obj.options.zonotopeOrder = zonotopeOrder;
            
        end
        
        % set intermediate order
        function set_intermediateOrder(obj, intermediateOrder)
            if intermediateOrder < 1
                error('Invalid intermediate order');
            end
            
            obj.options.intermediateOrder = intermediateOrder;
            
        end
        
        % set reduction technique
        function set_reductionTechnique(obj, reductionTechnique)
            obj.options.reductionTechnique = reductionTechnique;
        end
        
        % set error order
        function set_errorOrder(obj, errorOrder)
            if errorOrder < 1
                error('Invalid error order');
            end
            
            obj.options.errorOrder = errorOrder;
        end
        
        % set polytope Order
        function set_polytopeOrder(obj, polytopeOrder)
            
            if polytopeOrder < 1
                error('Invalid polytope order');
            end
            
            obj.options.polytopeOrder = polytopeOrder;
        end
        
        % set reduction Interval
        function set_reductionInterval(obj, reductionInterval)
            
            if reductionInterval < 0
                error('Invalid reduction interval');
            end
            
            obj.options.reductionInterval = reductionInterval;
            
        end
        
        % set max Error       
        function set_maxError(obj, maxError)
            
            if size(maxError, 1) ~= obj.dim && size(maxError, 2) ~= 1
                error('Invalid max error');
            end
            
            obj.options.maxError = maxError;
            
        end
        
        % set tensor Order
        function set_tensorOrder(obj, tensorOrder)
            if tensorOrder < 1
                error('Invalid tensor Order');
            end
            
            obj.options.tensorOrder = tensorOrder;
            
        end
        
        % set originContained
        function set_originContained(obj, originContained)
            obj.options.originContained = originContained;
        end
        
        % set advancedLinErrorComp
        function set_advancedLinErrorComp(obj, advancedLinErrorComp)
            obj.options.advancedLinErrorComp = advancedLinErrorComp;
        end
        
        % set plot Type
        function set_plotType(obj, plotType)
            obj.options.plotType = plotType;
        end
        
        % set uTrans
        function set_uTrans(obj, uTrans)
            obj.options.uTrans = uTrans;
        end
        
        % set input set U
        function set_U(obj, U)
            obj.options.U = U;
        end
        
        % set initial set
        function set_R0(obj, R0)
            obj.options.R0 = R0;
        end
        
        % set tFinal
        function set_tFinal(obj, tFinal)
            if tFinal <= 0
                error('tFinal should be > 0');
            end
            obj.options.tFinal = tFinal;
        end
        
        % set tStart
        function set_tStart(obj, tStart)
            
            if tStart < 0
                error('Invalid tStart');
            end
            
            obj.options.tStart = tStart;
        end
        
        % set inital state for simulation x0
        function set_x0(obj, x0)
            if length(x0) ~= obj.dim
                error('Dimension mismatch between initial state and the system');
            end
            
            obj.options.x0 = x0;
        end
        
    end
    
    
    % reachability anlaysis method
    methods
        
    end
    
end

