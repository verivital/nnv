classdef SignLayer < handle
    % Activation is a class that applies the given activation function to the
    % given input without changing its size
    % author: Mykhailo Ivashchenko
    % date: 09/052021
    
    properties
        f = 'sign'; % activation function;
        gamma = 0; % used only for leakyReLU layer
        mode = ""; % used for sign application mode
        
        option = []; % parallel option, 'parallel' or []
        dis_opt = []; % display option, 'display' or []
        lp_solver = 'linprog'; % lp solver option, 'linprog' or 'glpk'
        relaxFactor = 0; % use only for approx-star method
    end
    
    methods % constructor - evaluation - sampling
        % Constructor
        function obj = SignLayer(varargin)
            
            switch nargin
                case 0
                    obj.gamma = 0;
                    obj.mode = 'polar_zero_to_pos_one';
                case 1
                    obj.gamma = varargin{1};
                    obj.mode = 'polar_zero_to_pos_one';
                case 2
                    obj.gamma = varargin{1};
                    obj.mode = varargin{2};
                otherwise
                    error('Invalid number of input arguments, should be 1');
            end
            
            if obj.gamma >= 1
                error('Invalid parameter for leakyReLU, gamma should be <= 1');
            end
               
        end
        
        % Evaluation method
        function y = evaluate(obj, y1)  % evaluation of this layer with a specific vector    

            if ~strcmp(obj.f, 'sign')
                error('This class supports only the sign activation function.');
            end
            
            y = Sign.evaluate(y1, obj.mode);
        end
        
        % evaluate the value of the layer output with a set of vertices
        function Y = sample(obj, V)
  
            n = size(V, 2);
            Y = [];
            for j=1:n
                y = obj.evaluate(V(:, j));
                Y = [Y y];
            end
        end
        
    end
    
    
    methods % reachability analysis method
        
        function S = reach(varargin)
            % @I: an array of inputs
            % @method: 'exact-star' or 'approx-star' or 'approx-zono' or 'abs-dom', i.e., abstract domain (support
            % later) or 'face-latice', i.e., face latice (support later)
            % @option:  'parallel' use parallel computing
            %           '[]' or not declared -> don't use parallel
            %           computing
            
            % author: Dung Tran
            % date: 27/2/2019
            % update: 7/10/2020 add the relaxed approx-star method for poslin, logsig and tansig
            %         7/18/2020: add dis_play option + lp_solver option
             
            % parse inputs 
            switch nargin
                
                
                case 7
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    obj.option = varargin{4};
                    obj.relaxFactor = varargin{5}; % only use for approx-star method
                    obj.dis_opt = varargin{6};
                    obj.lp_solver = varargin{7};
                case 6
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    obj.option = varargin{4};
                    obj.relaxFactor = varargin{5}; % only use for approx-star method
                    obj.dis_opt = varargin{6};
                case 5
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    obj.option = varargin{4};
                    obj.relaxFactor = varargin{5}; % only use for approx-star method
                case 4
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    obj.option = varargin{4};
                case 3
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end
            
            if ~strcmp(method, 'exact-star') && ~strcmp(method, 'approx-star') && ~strcmp(method, 'approx-star-fast') && ~strcmp(method, 'approx-zono') && ~strcmp(method, 'abs-dom') && ~strcmp(method, 'exact-polyhedron') && ~strcmp(method, 'approx-star-split') && ~strcmp(method,'approx-star-no-split')
                error('Unknown reachability analysis method');
            end
            
            if strcmp(method, 'exact-star') && (~strcmp(obj.f, 'purelin') && ~strcmp(obj.f, 'leakyrelu') && ~strcmp(obj.f, 'poslin') && ~strcmp(obj.f, 'sign') && ~strcmp(obj.f, 'satlin') && ~strcmp(obj.f, 'satlins'))
                method = 'approx-star';
                fprintf('\nThe current layer has %s activation function -> cannot compute exact reachable set for the current layer, we use approx-star method instead', obj.f);
            end
            
            n = length(I);
            S = [];
%             gamma1 = obj.gamma;
            
            I1 = I;
            
            if obj.mode == ""
                obj.mode = 'polar_zero_to_pos_one';
            end
            
            if strcmp(obj.option, 'parallel') % reachability analysis using star set
                
                rF = obj.relaxFactor;
                dis = obj.dis_opt;
                lps = obj.lp_solver;
                mode_ = obj.mode;
                parfor i=1:n
                    S = [S Sign.reach(I1(i), method, [], rF, dis, lps, mode_)];
                end
            else
                
                for i=1:n
                    S = [S Sign.reach(I1(i), method, [], obj.relaxFactor, obj.dis_opt, obj.lp_solver, obj.mode)];
                end
            end
        end
    end    
end


