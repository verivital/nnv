classdef Operation
    % Operation : class for specific operation in Feedforward neural
    % network reachability
    
    % An Operation can be: 
    %   1) AffineMap
    %   2) PosLin_stepExactReach
    %   3) PosLin_approxReachZono
    %   4) PosLin_approxReachAbsDom
    %   5) PosLin_approxReachStar
    
    %   6) SatLin_approxReachStar
    %   7) SatLin_stepExactReach
    %   8) SatLin_approxReachZono
    %   9) SatLin_approxReachAbsDom
    
    %   10) SatLins_stepExactReach
    %   11) SatLins_approxReachStar
    %   12) SatLins_approxReachZono
    %   13) SatLins_approxReachAbsDom
    
    %   14) LogSig_approxReachZono
    %   15) LogSig_approxReachStar
    
    %   16) TanSig_approxReachZono
    %   17) TanSig_approxReachStar
    
    % The main method for Operation class is: 
    %   2) execute
    
    % The Operation class is used for verification of FFNN using Deep First
    % Search algorithm
    
    
    % Author: Dung Tran
    % Date: 1/18/2020
    
    properties
        
        Name = '';
        map_mat = []; % used if the operation is an affine mapping operation
        map_vec = []; % used if the operation is an affine mapping operation
        
        index = []; % used if the operation is PosLin or SatLin or SatLins stepReach operation
        % index is the index of the neuron the stepReach is performed
        
        method = ''; % reach method 
        
    end
    
    methods
        
        % constructor
        function obj = Operation(varargin)
            % @name: name of the operation
            % @W: affine mapping matrix 
            % @b: affine mapping vector
            % @index: neuron index
            
            % author: Dung Tran
            % date: 1/18/2020
            
            switch nargin
                
                case 3 
                    name = varargin{1};
                    mat = varargin{2};
                    vec = varargin{3}; 
                    
                    if ~strcmp(name, 'AffineMap')
                        error('The operation is not an affine mapping operation');
                    end
                    
                    if size(mat, 1) ~= size(vec, 1)
                        error('Inconsistent dimension between that affine mapping matrix and vector');
                    end
                    
                    if size(vec, 2) ~= 1
                        error('Affine mapping vector should have one column');
                    end
                    
                    obj.Name = name; 
                    obj.map_mat = mat;
                    obj.map_vec = vec; 
                
                case 2
                    
                    name = varargin{1};
                    index = varargin{2};
                    
                    if ~strcmp(name, 'PosLin_stepExactReach') && ~strcmp(name, 'SatLin_stepExactReach') && ~strcmp(name, 'SatLins_stepExactReach')   
                        error('Unknown operation name');
                    end
                    
                    if index < 1 
                        error('Invalid neuron index');
                    end
                    
                    obj.Name = name;
                    obj.index = index;
                    
                case 1 
                    
                    name = varargin{1};
                    
                    S1 = ~strcmp(name, 'PosLin_approxReachStar') && ~strcmp(name, 'PosLin_approxReachZono') && ~strcmp(name, 'PosLin_approxReachAbsDom');
                    S2 = ~strcmp(name, 'SatLin_approxReachStar') && ~strcmp(name, 'SatLin_approxReachZono') && ~strcmp(name, 'SatLin_approxReachAbsDom');
                    S3 = ~strcmp(name, 'SatLins_approxReachStar') && ~strcmp(name, 'SatLins_approxReachZono') && ~strcmp(name, 'SatLins_approxReachAbsDom');
                    S4 = ~strcmp(name, 'LogSig_approxReachStar')&& ~strcmp(name, 'LogSig_approxReachZono');
                    S5 = ~strcmp(name, 'TanSig_approxReachStar')&& ~strcmp(name, 'TanSig_approxReachZono');
                    
                    if  S1 && S2 && S3 && S4 && S5
                        error('Unknown operation name');
                    end
                    
                    obj.Name = name;
                    
                case 0
                    
                    
                otherwise
                    
                    error('Invalid number of arguments');
                    
            end
                         
        end
        
        
        % execute the operation
        function S = execute(obj, I)
            % @I: a star input set 
            % @S: a star output set or an array of star output sets 
            
            % author: Dung Tran
            % date: 1/18/2020
                        
                        
            if strcmp(obj.Name, 'AffineMap')
                S = I.affineMap(obj.map_mat, obj.map_vec);
            % PosLin
            elseif strcmp(obj.Name, 'PosLin_stepExactReach')
                [xmin, xmax] = I.estimateRange(obj.index); 
                S = PosLin.stepReach(I, obj.index, xmin, xmax);
            elseif strcmp(obj.Name, 'PosLin_approxReachStar')
                S = PosLin.reach_star_approx2(I);
            elseif strcmp(obj.Name, 'PosLin_approxReachZono')
                S = PosLin.reach_zono_approx(I);
            elseif strcmp(obj.Name, 'PosLin_approxReachAbsDom')
                S = PosLin.reach_abstract_domain(I);
            % SatLin
            elseif strcmp(obj.Name, 'SatLin_stepExactReach')
                S = SatLin.stepReach(I, obj.index);
            elseif strcmp(obj.Name, 'SatLin_approxReachStar')
                S = SatLin.reach_star_approx(I);
            elseif strcmp(obj.Name, 'SatLin_approxReachZono')
                S = SatLin.reach_zono_approx(I);
            elseif strcmp(obj.Name, 'SatLin_approxReachAbsDom')
                S = SatLin.reach_abstract_domain(I);
            % SatLins
            elseif strcmp(obj.Name, 'SatLins_stepExactReach')
                S = SatLins.stepReach(I, obj.index);
            elseif strcmp(obj.Name, 'SatLins_approxReachStar')
                S = SatLins.reach_star_approx(I);
            elseif strcmp(obj.Name, 'SatLin_approxReachZono')
                S = SatLins.reach_zono_approx(I);
            elseif strcmp(obj.Name, 'SatLins_approxReachAbsDom')
                S = SatLins.reach_abstract_domain(I);
            % LogSig    
            elseif strcmp(obj.Name, 'LogSig_approxReachStar')
                S = LogSig.reach_star_approx(I);
            elseif strcmp(obj.Name, 'LogSig_approxReachZono')
                S = LogSig.reach_zono_approx(I);
            % TanSig
            elseif strcmp(obj.Name, 'TanSig_approxReachStar')
                S = TanSig.reach_star_approx(I);
            elseif strcmp(obj.Name, 'TanSig_approxReachZono')
                S = TanSig.reach_zono_approx(I);
            else
                error('Unknown operation')
            end
            
            
        end
        
        
        
    end
    
end

