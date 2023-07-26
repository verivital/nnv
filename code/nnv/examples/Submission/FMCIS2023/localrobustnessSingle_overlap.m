function robust_overlap = localrobustnessSingle_overlap(LB,UB,allowable_LB,allowable_UB)
%% the function 'localrobustnessSingle_overlap' takes the following inputs
%   LB              : lower bound from the output reachable set
%   UB              : upper bound from the output reachable set
%   allowbale_LB    : permissible lower bound
%   allowable_UB    : permissible upper bound
% and outputs the Percentage Overlap i.e., the overlap between the 
% estimated range and the permissible range for a single time-instance

    if  UB <= allowable_UB && LB >= allowable_LB
        robust_overlap = 1;
        %fprintf('case 1');
    elseif  UB < allowable_LB || LB > allowable_UB
        robust_overlap = 0;
    %     fprintf('case 2');
    elseif UB < allowable_UB && UB > allowable_LB && LB < allowable_LB
        robust_overlap = (UB - allowable_LB)/(UB - LB);
    %     fprintf('case 3a');
    elseif UB > allowable_UB && LB > allowable_LB && LB < allowable_UB
        robust_overlap = (allowable_UB - LB)/(UB - LB);
    %     fprintf('case 3b');
    else
        robust_overlap = (allowable_UB - allowable_LB)/(UB - LB);
    %     fprintf('case 4'); 
    end
end