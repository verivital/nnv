classdef NeuralODE < handle
    % NeuralODE Class is a class for Verification of Neural Networks with
    % NeuralODE layers
    % Author: Diego Manzanas
    % Date: 01/18/2019
    
    properties
        
        Name = 'nnode'; % name of the network
        Layers = {}; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        numLayers = 0; % number of Layers
%         numNeurons = 0; % number of Neurons
        InputSize = 0; % number of Inputs
        OutputSize = 0; % number of Outputs
        
        % properties for reach set computation
        
        reachMethod = 'approx-star';    % reachable set computation scheme, default - 'approx-star'
        relaxFactor = 0; % default - solve 100% LP optimization for finding bounds in 'approx-star' method
        reachOption = []; % parallel option, default - non-parallel computing
        numCores = 1; % number of cores (workers) using in computation
        reachSet = [];  % reachable set for each layers
        outputSet = [];
        reachTime = []; % computation time for each layers
        totalReachTime = 0; % total computation time
        
        features = {}; % outputs of each layer in an evaluation
        dis_opt = []; % display option = 'display' or []
        lp_solver = 'glpk'; % choose glpk as default LP solver for constructing reachable set
        node_type = 'classification'; % Two types of Neural ODE supported: 'classification', 'timeseries'
        
    end
    
    
    methods % constructor, evaluation, sampling, print methods
        
        % constructor
        function obj = NeuralODE(varargin)
            
            switch nargin
                case 4
                    name = varargin{1};
                    Layers = varargin{2};
                    inputsize = varargin{3};
                    outputsize = varargin{4};
                    nL = length(Layers); % number of Layers
                                        
                    obj.Name = name;
                    obj.Layers = Layers;
                    obj.numLayers = nL;    % number of layers
                    obj.InputSize = inputsize; % input size
                    obj.OutputSize = outputsize; % output size
                    
                case 1
                    obj.Name = 'NeuralODE';
                    obj.Layers = varargin{1};
                    obj.numLayers = length(obj.Layers);
                    obj.InputSize = 0;
                    obj.OutputSize = 0;
                    
                case 0
                    
                    obj.Layers = {};
                    obj.numLayers = 0;
                    obj.InputSize = 0;
                    obj.OutputSize = 0;
                    
                otherwise
                    error('Invalid number of inputs, should be 0,1 4');
            end
            
            % Check if neural ode is timeseries or classification
            for l = obj.Layers
                if contains(class(l{1}),'ODE')
                    if l{1}.time_series
                        obj.node_type = 'timeseries';
                    end
                end
            end
            
            % Compute number of inputs and outputs of Neural ODE
            
        end
                
        
        % Evaluation of a NeuralODE
        function y = evaluate(obj, x)
            % Evaluation of this NeuralODE
            % @x: input vector x
            % @y: output vector y
            % @features: output of all layers
            
%             y = x;
            for i=1:obj.numLayers
                xt = [];
                for yp = 1:size(x,2)
                    y = obj.Layers{i}.evaluate(x(:,yp));
                    xt = [xt y];
                end
                obj.features{i} = xt;
                x = xt;
            end
            y = x;
        end
        
    end
    
    
    methods % reachability analysis method
        
        function [IS, reachTime] = reach(varargin)
            % [IS, reachTime] = reach(varargin)
            % IS = reach output set
            % reachTime = total computation time
            % INPUTS
            % @I: input set, a star set          
            % @method: = 'exact-star' or 'approx-star' -> compute reach set using stars
            
            switch nargin 
                
                case 2
                    
                    obj = varargin{1};
                    if ~isstruct(varargin{2})
                        inputSet = varargin{2};
                    else
                       if isfield(varargin{2}, 'inputSet')
                           inputSet = varargin{2}.inputSet;
                       else
                           error('No input set for reachability analysis');
                       end
                       if isfield(varargin{2}, 'reachMethod')
                           obj.reachMethod = varargin{2}.reachMethod;
                       end
                       if isfield(varargin{2}, 'numCores')
                           obj.numCores = varargin{2}.numCores;
                       end
                       if isfield(varargin{2}, 'relaxFactor')
                           obj.relaxFactor = varargin{2}.relaxFactor;
                       end
                       if isfield(varargin{2}, 'dis_opt')
                           obj.dis_opt = varargin{2}.dis_opt; % use for debuging
                       end
                    end
                    
                case 3 
                    
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = 1; 
                    
                case 4
                    
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = varargin{4};
                    
                case 5
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = varargin{4};
                    obj.relaxFactor = varargin{5};
                    
                case 6
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = varargin{4};
                    obj.relaxFactor = varargin{5};
                    obj.dis_opt = varargin{6}; % use for debuging
                    
                otherwise 
                    
                    error('Invalid number of input arguments, the number should be 1, 2, 3, 4, 5');
                
            end
                       
            if  obj.numCores > 1 % Parallel computation not supported for neural ODEs
%                 obj.start_pool;
%                 obj.reachOption = 'parallel';
                obj.reachOption = [];
                disp('Parallel computation not supported for Neural ODEs');
            else
                obj.reachOption = [];
            end
            
            obj.reachSet = cell(1, obj.numLayers);
            obj.reachTime = zeros(1, obj.numLayers);
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for the network %s...', obj.Name);
            end
            rs = inputSet;
            for i=2:obj.numLayers+1
                if strcmp(obj.dis_opt, 'display')
                    fprintf('\nPerforming analysis for Layer %d (%s)...', i-1, obj.Layers{i-1}.Name);
                end
                start_time = tic;
                rs_new = obj.Layers{i-1}.reach(rs, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, 'glpk');
                obj.reachTime(i-1) = toc(start_time);
                rs = rs_new;
                obj.reachSet{i-1} = rs_new;
                if strcmp(obj.dis_opt, 'display')
                    fprintf('\nReachability analysis for Layer %d (%s) is done in %.5f seconds', i-1, obj.Layers{i-1}.Name, obj.reachTime(i-1));
                    fprintf('\nThe number of reachable sets at Layer %d (%s) is: %d', i-1, obj.Layers{i-1}.Name, length(rs_new));
                end
            end
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nReachability analysis for the network %s is done in %.5f seconds', obj.Name, sum(obj.reachTime));
                fprintf('\nThe number ImageStar in the output sets is: %d', length(rs_new));
            end
            obj.totalReachTime = sum(obj.reachTime);
            IS = rs_new;
            obj.outputSet = rs_new;
            reachTime = obj.totalReachTime;
        end
        
        
    end
    
    
    methods(Static)
        
        % check robustness using the outputSet
        function [rb, cands] = checkRobust(outputSet, correct_id)
            % @outputSet: the outputSet we need to check
            % @correct_id: the correct_id of the classified output
            % @rb: = 1 -> robust
            %      = 0 -> not robust
            %      = 2 -> unknown
            % @cand: possible candidates
            
            R = outputSet;
            
            if contains(class(outputSet),'Image')
                if correct_id > outputSet.numChannel || correct_id < 1
                    error('Invalid correct id');
                end

                R = R.toStar;
%                 [lb, ub] = R.estimateRanges;
            end
            [lb, ub] = R.getRanges;
            [~, max_ub_id] = max(ub);
            cands = [];
            if max_ub_id ~= correct_id
                rb = 2;
                cands = max_ub_id;
            else                   
                
                max_val = lb(correct_id);
                max_cd = find(ub > max_val); % max point candidates
                max_cd(max_cd == correct_id) = []; % delete the max_id

                if isempty(max_cd)
                    rb = 1;
                else            
                    
                    n = length(max_cd);
                    C1 = R.V(max_cd, 2:R.nVar+1) - ones(n,1)*R.V(correct_id, 2:R.nVar+1);
                    d1 = -R.V(max_cd, 1) + ones(n,1)*R.V(correct_id,1);
                    S = Star(R.V, [R.C;C1], [R.d;d1], R.predicate_lb, R.predicate_ub);
                    if S.isEmptySet
                        rb = 2;
                        cands = max_cd;
                    else                       
                        count = 0;
                        for i=1:n
                            if R.is_p1_larger_than_p2(max_cd(i), correct_id)
                                rb = 2;
                                cands = max_cd(i);
                                break;
                            else
                                count = count + 1;
                            end
                        end
                        if count == n
                            rb = 1;
                        end         
                    end    

                end

            end

        end
        
        
        
    end
    
    
    
end

