function out = join(obj)
%
% merges unions of polyhedra in the same dimension
%

% kick out empty polyunions
dm = cell(size(obj));
[dm{:}] = obj.Dim;
obj(cellfun('isempty',dm)) = [];
obj([dm{:}]==0);

no = numel(obj);
if no<=1
    % nothing to do
    out = obj;
else    
    dims = zeros(no, 1);
    for i=1:no
        dims(i) = obj(i).Dim;
    end
    
    if any(diff(dims))
        error('Only polyunions of the identical dimensions can be joined.');
    end

    % do not check convexity, overlaps, and connectivity
    % check boundedness
    bnd = cell(no, 1);
    for i=1:no
        bnd{i} = obj(i).Internal.Bounded;
    end
    ie = cellfun('isempty',bnd);
    if any(~ie)
        % one of unions has Bounded property set
        if any([bnd{:}]==1)
            % check for bounded sets
            b = obj.isBounded;
            if ~all(b)
                error('All unions must be bounded.')
            end
        end
        indexb = find(~ie);
        isb = bnd{indexb(1)};
    else
        isb = [];
    end
    
    % check dimensionality
    d = cell(no, 1);
    for i=1:no
        d{i} = obj(i).Internal.FullDim;
    end
    id = cellfun('isempty',d);
    if any(~id)
        % one of unions has FullDim property set
        if any([d{:}]==1)
            % check for fulldimensionality
            df = obj.isFullDim;
            if ~all(df)
                error('All unions must be full-dimensional.')
            end
        end
        indexfd = find(~id);
        isfd = d{indexfd(1)};
    else
        isfd = [];
    end
    
    out = PolyUnion('Set', obj(1).Set);
    if numel(obj)>1
        for i = 2:length(obj)
            out.add(obj(i).Set);
        end
    end
    
    % set internal properties
    out.Internal.Bounded = isb;
    out.Internal.FullDim = isfd;
    
end


end
