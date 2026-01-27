classdef CSV
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        
        % get the conservativeness
        
        function [R, r] = getConservativeness(approx_range, actual_range)
            % this function returns the conservativeness of
            % over-approximate reachable set
            % @approx_range: the approximate range of the reachable set
            %                contain lower-bound and upper-bound vector of
            %                the output set
            % @actual_range: the range of the actual reachable set
            %                = [lb ub] the lower-bound and upper-bound
            %                vector of the actual reachable set
            
            % @R: conservativeness in %
            % @r: = 0 : exact reachable set
            %     = 1 : the approximate reach set is an over-approximation of the actual reachable set
            %     = 2 : the approximate reach set is an under-approximation of the actual reachable set
            %     = 3 : the approximate reachable set is wrong, it is
            %     neigther an over-approximation or under-approximation of
            %     the actual reachable set.
            
            [n1, m1] = size(approx_range);
            [n2, m2] = size(actual_range);
            
            if n1 ~= n2 || m1 ~= m2
                error('Inconsitent dimension between approximate range and actual range');
            end
            
            if m1 ~= 2 || m2 ~= 2
                error('approx_range or actual range should be n x 2 matrix');
            end
            
            
            % compute the conservativeness
            R = zeros(n1, 1);
            for i=1:n1
                R(i) = 100 * max(abs(approx_range(i, 1) - actual_range(i,1)), abs(approx_range(i,2)-actual_range(i,2)))/(actual_range(i,2) - actual_range(i,1));
            end
  
            r = 0;
            temp1 = 0;
            temp2 = 0;

            for i=1:n1
                if R(i) >= 0.01                    
                    if approx_range(i,1) <= actual_range(i,1) && approx_range(i, 2) >= actual_range(i, 2)
                      temp1 = temp1 +1;
                    end
                    
                    if approx_range(i, 1) >= actual_range(i,1) && approx_range(i, 2) <= actual_range(i, 2)
                        temp2 = temp2 + 1;
                    end
   
                end
            end
            
            if temp1 == n1
                r = 1; % approximate set is an over-approximation of the actual reachable set
            elseif temp2 == n1
                r = 2; % approximate set is an under-approximation of the actual reachable set
            elseif sum(R) >= 0.01
                r = 3;
            end
            
            
           
            
        end
        
    end
    
end

