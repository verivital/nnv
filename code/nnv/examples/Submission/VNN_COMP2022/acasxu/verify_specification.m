function [result] = verify_specification(reachSet, property)
%VERIFY_SPEFICIATION verify a vnnlib-defined property (property) based on the output
%ranges of a neural network (reachSet)
% Syntax:
%    [result] = verify_speficication(reachSet, property)
% Inputs:
%   - reachSet: computed output set of neural network (e.g. 1x1 Star)
%   - property: cell array of cell array defining all conditions that outputSet must satisfy
% Output:
%   - result: 0 ->  property failed
%             1 ->  property satisfied
%             2 ->  unknown

[r_lb, r_ub] = reachSet.getRanges(); % get ranges of reach set
i = 1;
condition = string();
while i <= length(properties) % we are going to use eval('expression') to evaluate the properties
    % Idea: 
    % convert entire property to string, then run it inside eval. May need to do this 
    % iteratively to add the || and && conditions, that may be a little more complex.
    
end % end while loop

end % close function

