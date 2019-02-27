classdef SatLin
    % SATLIN : class for computing reachable set of Satlin Transfer Function 
    %   Reference: https://www.mathworks.com/help/deeplearning/ref/satlin.html
    % Author: Dung Tran
    % Date: 27/2/2019
    
    properties
        
    end
    
    methods(Static)
        
        % stepSatLin method, compute reachable set for a single step
        function S = stepReach(I, index)
            % @I: single star set input
            % @index: index of the neural performing stepSatLin
            % @S: star output set
            
            % author: Dung Tran
            % date: 27/2/2019
            
            
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
        
        
        
        
    end
end

