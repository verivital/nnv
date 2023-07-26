function POR = globalrobustnessSingle_overlap(LB,UB,allowable_LB,allowable_UB)
%% the function 'globalrobustnessSingle_overlap' takes the following inputs
%   LB              : lower bound from the output reachable set
%   UB              : upper bound from the output reachable set
%   allowbale_LB    : permissible lower bound
%   allowable_UB    : permissible upper bound
% and outputs the percentage overlap robustness measure as POR

    cR = 0;
    total = length(LB);
    for i = 1: total
        robust = localrobustnessSingle_overlap(LB(i),UB(i),allowable_LB(i),allowable_UB(i));
        cR = cR + robust;
    end
    POR = cR/total;
end