function PR = globalrobustnessSingle(LB,UB,allowbale_LB,allowable_UB)
%% the function 'globalrobustnessSingle' takes the following inputs
%   LB              : lower bound from the output reachable set
%   UB              : upper bound from the output reachable set
%   allowbale_LB    : permissible lower bound
%   allowable_UB    : permissible upper bound
% and outputs the percentage overlap robustness measure as PR

    cR = 0;
    total = length(LB);
    for i = 1: total
        robust = localrobustnessSingle(LB(i),UB(i),allowbale_LB(i),allowable_UB(i));
        cR = cR + robust;
    end
    PR = cR/total;
end