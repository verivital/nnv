function obj = add(obj, C)
%
%  Add the object to the union 
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);

if ~isa(C,'Polyhedron')
    error('Argument must be a Polyhedron object.');
end

% remove empty sets
c = isEmptySet(C);
C(c) = [];

if numel(C)<1
    return;
end

% deal with arrays
if numel(obj)>1
    for i=1:numel(obj)
        obj(i) = obj(i).add(C);
    end
    return;
end

% check dimension
D = [C.Dim];
if any(D(1)~=D)
    error('All polyhedra must be in the same dimension.');
end
Dim = D(1);
if obj.Dim~=Dim
    error('The polyhedra must be in the dimension %d.',obj.Dim);
end


% length of the potential elements
N = numel(C);

% reserve -1 statements to indicate that these properties have not been
% checked
% t = false means that the property is ok
% t = true means that the property is violated
s = {'Convex','Overlaps','Connected','Bounded','Fulldim'};
t = -ones(length(s),1);

% check convexity 
if ~isempty(obj.Internal.Convex)
    t(1) = false;
    if obj.Internal.Convex        
        % compute the convex hull
		Pn = obj.convexHull();
        
        % merge existing polyhedra Pn and new C
        S = [Pn(:); C(:)];
        
		% detect convexity
		t(1) = ~PolyUnion(S).isConvex();
    end
end

% check overlaps
% here the value of t(2) is opposite
if ~isempty(obj.Internal.Overlaps)
    t(2) = false;
    if ~obj.Internal.Overlaps 
        % if there exist convex hull and the union is convex, test intersection with this
        if t(1)==true && isfield(obj.Internal,'convexHull')            
            for i=1:N
                IH = intersect(obj.Internal.convexHull,C(i));
                if isFullDim(IH)
                    t(2) = true;
                    break;
                end
            end
        else
            % run a test consecutively for each polyhedron in the set
            for i=1:N
                IC = intersect(obj.Set,C(i));
                if any(IC.isFullDim)
                    t(2) = true;
                    break;
                end
            end

        end
    end    
end

% check connectivity
if ~isempty(obj.Internal.Connected)
    t(3) = false;
    if obj.Internal.Connected
        
        % if the test for convexity passed, we don't need to check for
        % connectivity
        if t(1)~=false

            % test if the both sets contain one common point
            for i=1:N
                IC = intersect(obj.Set,C(i));
                for j=1:numel(IC)
                    xc = IC(j).interiorPoint;
                    if ~isempty(xc.x)                        
                        % is connected, stop search
                        bk = false;
                        break
                    else
                        bk = true;
                    end
                end
                if bk
                    % stop search
                    break;
                end
            end
            t(3)=bk;
        end
    end
end
    
% check boundedness
if ~isempty(obj.Internal.Bounded)
    t(4) = false;
    if obj.Internal.Bounded
        b = C.isBounded;
        if ~all(b)
            t(4) = true;
        end
    end    
end

% check full dimensionality 
if ~isempty(obj.Internal.FullDim)
    t(5) = false;
    if obj.Internal.FullDim
        f = C.isFullDim;
        if ~all(f)
            t(5) = true;
        end
    end
end


% if any statement is violated, throw an error
if any(t==true)
    str = '';
    iz = find(t==true)';
    for i=iz
        if i==iz(1)
            str = [str,s{i}];
        else
            str = [str,', ',s{i}];
        end
    end
    error('The polyhedra cannot be added because it conflicts with "%s" property.',str);
else
    % check functions
	
	fnames = obj.listFunctions;
	nf = length(fnames);
	for i = 1:N
		if any(~C(i).hasFunction(fnames))
			error('The set %i to be added holds different function names than the union.',i);
		elseif length(C(i).Functions) ~= nf
			error('All sets to be added must have associated %d number of functions.',nf);
		end
	end
    
    % if all properties are satisfied, it is ok to add the set to the union
    if obj.Num==0
        obj.Set = C(:);
    else
        obj.Set(obj.Num+1:obj.Num+N,1) = C(:);
    end
    
    % update the convex hull if it was computed
    if exist('H','var')
        % update the convexhull info
        obj.Internal.convexHull = H;
    end

end
      

end
