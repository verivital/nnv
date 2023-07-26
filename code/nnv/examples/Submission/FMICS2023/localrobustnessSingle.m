function robust = localrobustnessSingle(LB,UB,allowbale_LB,allowable_UB)
%% the function 'localrobustnessSingle' takes the following inputs
%   LB              : lower bound from the output reachable set
%   UB              : upper bound from the output reachable set
%   allowbale_LB    : permissible lower bound
%   allowable_UB    : permissible upper bound
% and outputs the robustness value i.e., if the reachable bounds fall within
% the permissible range or not, for a single time-instance

    if LB > allowbale_LB && UB < allowable_UB
        robust = 1;
    else 
        robust = 0;
    end
end