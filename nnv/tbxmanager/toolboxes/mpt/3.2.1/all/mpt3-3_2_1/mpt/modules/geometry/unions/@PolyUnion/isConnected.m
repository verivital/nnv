function ts = isConnected(obj)
%
% check if the polyhedron array forms a connected union
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% deal with arrays
if numel(obj)>1
    ts = -ones(size(obj));
    for i=1:numel(obj)
        ts(i) = obj(i).isConnected;
    end
    return;
end

% empty obj
if obj.Num==0
    ts = false;
    return;
elseif obj.Num==1
	% single region is connected
	ts = true;
	return
end

% check connectivity
if isempty(obj.Internal.Connected)        
    % if there is info about the convexity, use this info to determine
    % connectivity
    if ~isempty(obj.Internal.Convex)        
        if obj.Internal.Convex
            % if the union is convex, then it is connected as well
            % otherwise we have to check for connectivity
            ts = true;
            obj.Internal.Connected = ts;
            return
        end
    end
    
    % Each polyhedron inside the array must contain a common point
    % with any of the remaining sets.
    % The easiest would be to check connectivity between neighbors,
    % but since we don't know what are the neighbors, we run a
    % consecutive test for each polyhedron in the array reducing of
    % those polyhedra that have been checked.
    
    % get all combinations with
    c = cbn(1:obj.Num);
    
    % run a test consecutively for each combination of polyhedra in the set
    t = false(1,obj.Num);
    while size(c,1)>0
        if MPTOPTIONS.verbose >= 1
            fprintf('Regions: %d-%d \n',c(1,1),c(1,2));
        end
        
        % take always the first combination
        reg = c(1,:);
        IC = intersect(obj.Set(reg(1)),obj.Set(reg(2)));
        if ~IC.isEmptySet
            % region c(i,1)=reg(1)  is connected with c(i,2)=reg(2)
            t(reg(1)) = true;
            t(reg(2)) = true;
            % kick out region c(i,1) and c(i,2) from the list
            if ~isempty(c)
                c(c(:,1)==reg(1),:) = [];
            end
            if ~isempty(c)
                c(c(:,1)==reg(2),:) = [];
            end
        else
            % this combination has been checked, remove it from the list
            c(1,:)=[];
        end
        
    end
    % if all are connected, return true
    ts = all(t);
    obj.Internal.Connected = ts;
    
else
   ts = obj.Internal.Connected; 
end



end

function y = cbn(v)
%
% create two combinations of sets to be checked
% v - vector of indices
%

N = numel(v);

y = zeros(N*(N-1),2);
k=1;
for i=1:N
    for j=setdiff(1:N,i)
        y(k,:) = [v(i),v(j)];
        k = k+1;
    end
end



end
