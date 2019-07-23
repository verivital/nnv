function ts = isOverlapping(obj)
%
% check if the polyhedron array is overlapping or not
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% deal with arrays
if numel(obj)>1
    ts = -ones(size(obj));
    for i=1:numel(obj)
        ts(i) = obj(i).isOverlapping;
    end
    return;
end

% empty obj
if obj.Num <= 1
	% no regions = no overlaps
	% single region = no overlaps
    ts = false;
    return;
end

if isempty(obj.Internal.Overlaps)
    ts = false;
    % get all combinations without repetitions
    % the function "combnk" comes from statistical toolbox
    %c = combnk(1:obj.Num,2);
    c = nchoosek(1:obj.Num,2);
    
    % run a test consecutively for each combination of polyhedra in the set    
    for i=1:size(c,1)
        if MPTOPTIONS.verbose >= 1
            fprintf('Regions: %d-%d \n',c(i,1),c(i,2));
        end

        IC = intersect(obj.Set(c(i,1)),obj.Set(c(i,2)));
        if ~obj.Set(c(i,1)).isFullDim || ~obj.Set(c(i,2)).isFullDim
            if ~IC.isEmptySet
                ts = true;
                break;
            end
        else
            if IC.isFullDim
                ts = true;
                break;
            end
        end
    end
    obj.Internal.Overlaps = ts;
else
    ts = obj.Internal.Overlaps;
end



end
