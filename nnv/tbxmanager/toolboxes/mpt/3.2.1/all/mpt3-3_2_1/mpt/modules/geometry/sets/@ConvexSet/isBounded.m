function ts = isBounded(obj)
%
% Tests if the set is bounded.
% <p>
% Algorithm:
% <ol>
% <li>Compute the support in the positive and negative directions of each elementary vector.</li>
% <li>If any are unbounded, then the set is unbounded, otherwise not. </li>
% </ol>
%
% @return <code>true</code> if bounded, <code>false</code> otherwise
% true if empty
%

no = numel(obj);
if no>1
    ts = true(size(obj));
    for i=1:no
        ts(i) = obj(i).isBounded;
    end
    return
end

% empty set -> bounded
if obj.isEmptySet,
    ts = true;
    return;
end

% Compute the support in all +- elementary directions.
% Bounded iff all are bounded
ts = true;
I = eye(obj.Dim);
for i=1:obj.Dim
    if isinf(obj.support( I(i,:)'))
        ts = false;
        return;
    end
    if isinf(obj.support(-I(i,:)')),
        ts = false;
        return;
    end
end

end
