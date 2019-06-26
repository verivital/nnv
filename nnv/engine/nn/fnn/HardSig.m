classdef HardSig
    % HardSig : class for computing reachable set of HardSigmoid Transfer Function 
    % Reference https://www.tensorflow.org/api_docs/python/tf/keras/backend/hard_sigmoid
    % Author: Dung Tran
    % Date: 4/6/2019
    
    properties

    end
    
    methods(Static)
        % evaluation
        function y = evaluate(x)
            n = length(x);
            y = zeros(n,1);
            for i=1:n
               if x(i) < 2.5
                   y(i) = 0;
               end
               if x(i) > 2.5
                   y(i) = 1;
               end
               if x >= -2.5 && x <= 2.5
                   y(i) = 0.2 * x(i) + 0.5;
               end
            end 
        end
            
        % stepReach method, compute reachable set for a single step
        function S = stepReach(I, index)
            % @I: single star set input
            % @index: index of the neural performing stepSatLin
            % @S: star output set
            
            % author: Dung Tran
            % date: 4/6/2019
            
            
            if ~isa(I, 'Star')
                error('Input is not a star set');
            end
            
            % case 1: x(index) <= 0, SatLin(x(index)) = 0
            C0 = I.V(index, 2:I.nVar + 1);
            d0 = I.V(index, 1);
            C1 = [I.C; C0];
            d1 = [I.d; -d0];
            V1 = I.V;
            V1(index, :) = zeros(1, I.nVar + 1);
            S1 = Star(V1, C1, d1);
            if S1.isEmptySet
                S1 = [];
            end
            
            % case 2: 0 <= x(index) <= 1, SatLin(x(index)) = x(index)            
            V2 = I.V;
            C2 = [I.C; C0; -C0];
            d2 = [I.d; 1-d0; d0];
            S2 = Star(V2, C2, d2);
            if S2.isEmptySet
                S2 = [];
            end
            
            % case 3: x(index) >= 1
            C3 = [I.C; -C0];
            d3 = [I.d; d0 - 1];
            V3 = I.V;
            V3(index, 1) = 1;
            V3(index, 2:I.nVar + 1) = zeros(1, I.nVar);
            S3 = Star(V3, C3, d3);
            if S3.isEmptySet
                S3 = [];
            end
            
            S = [S1 S2 S3];
                     
        end
        
        
        % stepReach with multiple inputs
        function S = stepReachMultipleInputs(varargin)
            % @I: an array of stars
            % @index: index where stepReach is performed
            % @option: = 'parallel' use parallel computing
            %          = not declare -> don't use parallel computing
            
            % author: Dung Tran
            % date: 27/2/2019
            
            switch nargin
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    option = varargin{3};
                case 2
                    I = varargin{1};
                    index = varargin{2};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2 or 3)');
            end
            
            
            
            p = length(I);
            S = [];
            
            if isempty(option)
                
                for i=1:p
                    S =[S, SatLin.stepReach(I(i), index)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    S =[S, SatLin.stepReach(I(i), index)];
                end
                
            else
                error('Unknown option');
            end
            
            
        end
       
        
        
        % exact reachability analysis using Star
        function S = reach_star_exact(I, option)
            % @I: an array of star input sets
            % @option: = 'parallel' use parallel option
            %          = '' do use parallel option
            
            % author: Dung Tran
            % date: 27/2/2019
            
            if ~isempty(I)       
                dim = I(1).dim;
                In = I;
                for i=1:dim
                    fprintf('\nPerforming SatLin_%d operation', i);
                    In = SatLin.stepReachMultipleInputs(In, i, option);
                end             
                
                S = In;
            else
                S = [];
            end
            
              
        end
        
    end
end

