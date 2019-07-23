function ts = isEmptySet(obj)
% ISEMPTY Tests if the set is empty.
%
% tf = isempty()
%
% Returns:
%   true if this set is empty, false otherwise
%

global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

no = numel(obj);
if no>1
    ts = false(size(obj));
    for i=1:no
        ts(i) = obj(i).isEmptySet;
    end
    return
elseif no<1
    ts = true;
    return
end

% Try to compute support - empty if this is infeasible
ret = obj.extreme(ones(obj.Dim,1));
ts = false;
if ret.exitflag == MPTOPTIONS.INFEASIBLE
    ts = true; 
end
end
