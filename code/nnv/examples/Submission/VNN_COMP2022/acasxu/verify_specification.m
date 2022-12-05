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

    [lb, ub] = reachSet.getRanges(); % get ranges of reach set
    disp(' ');
    i = 1;
    result = ones(length(property),1);
    for i = 1:length(property)
        result(i) = eval(property{i}{1}); % Verification result (1: sat, 0: unsat or unknown)
        if ~result(i) % check it property is sat
            result(i) = eval(property{i}{2}); % if result = 1, this means we prove unsat
            if ~ result(i)
                result(i) = 2; % unknown
            else
                result(i) = 0; % unsat
                break;
            end
        end
    end % end while loop
    % check the global result of the specifications
    if all(result == 1)
        result = 1;
    elseif any(result== 0)
        result = 0;
    else
        result = 2;
    end
end % close function

