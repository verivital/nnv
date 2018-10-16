classdef Verifier
    % Verifier class 
    %   Contains methods for checking safety property of a Neural network
    
    properties
    end
    
    methods(Static)
        
        % check safety
        function [safe, t] = checkSafety(reachSet, property)
            % @reachSet: the computed reachable set, is a set of polyhedra
            % @property: a set of safety properties, is a set of polyhedra
            
            t1 = tic;
            n = length(reachSet);
            m = length(property);
            
            for i=1:n
                if ~isa(reachSet(i), 'Polyhedron')
                    error('the %d^th reachable set is not a polyhedron', i);
                end
            end
            
            for j=1:m
                if ~isa(property(j), 'Polyhedron')
                    error('the %d^th property is not a polyhedron', j);
                end
            end
            
            N = 0;
            for i=1:n
                for j=1:m
                   % fprintf('\nchecking for i=%d, j=%d', i, j);
                   if ~isEmptySet(intersect(reachSet(i), property(j)))
                        fprintf('\nThe %d^th reachable set violates the %d^th property', i, j);
                        N = N + 1;
                    end
                end
            end
            
            if N == 0
                safe = true;
            else
                safe = false;
            end
            
            t = toc(t1);
            
        end
        
    end
    
end

