function approx = outerApprox(obj)
% OUTERAPPROX Computes an outer approximation of this set.
%
% -------------------------------------------------------------------
% approx = outerApprox() : Compute the smallest axis-aligned hypercube
%                          that contains this set.
%
% Returns:
%  approx - Polyhedron (Hypercube)
%  approx.Internal.lb / ub contain vectors of upper and lower bounds.
%
% -------------------------------------------------------------------
% approx = outerApprox(varargin) : Compute a containing polyhedron of
%                                  specified complexity or error.
%
% Implements the <a
% href="http://control.ee.ethz.ch/~cjones/geometricComputing.php">implicit double description method </a>
%
% Returns:
%  approx - Polyhedron
%
%
%
% @param maxFacets [default <eq>\infty</eq>] Maximum number of facets desired in the approximation
% @param maxVertices [default <eq>\infty</eq>] Maximum number of vertices desired in the approximation
% @param maxErr [default <eq>0</eq>] Maximum error in the Hausdorff metric between the approximation and the set
% @return Polytope containing the set of given complexity
%
% Input case 2:
%
%
% <p>
% Algorithm:<br>
% <code> for each i=1..d  support(-e_i) <= x_i <= support(e_i) </code><br>
% where e_i is the i^th unit vector
%
% @return Smallest axis-aligned hypercubes that contains this set.
%  or null if set is unbounded or empty
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if obj.isEmptySet
    approx = Polyhedron.emptySet(obj.Dim);
    approx.Internal.lb =  inf*ones(obj.Dim,1);
    approx.Internal.ub = -inf*ones(obj.Dim,1);
    return
end

% deal with arrays
no = numel(obj);
if no>1
    approx(size(obj)) = Polyhedron;
    for i=1:no
        approx(i) = obj(i).outerApprox;
    end
    return
end


if nargin < 2
    % Compute the support in all +- elementary directions.
    I   = eye(obj.Dim);
    pos = zeros(obj.Dim,1);
    neg = zeros(obj.Dim,1);
    for i=1:obj.Dim
        pos(i) =  obj.support( I(i,:)');
        neg(i) = -obj.support(-I(i,:)');
    end
    
    H = [eye(obj.Dim), pos; -eye(obj.Dim), -neg];
    
    approx = Polyhedron(H(:, 1:end-1), H(:, end));
    approx.Internal.lb = neg;
    approx.Internal.ub = pos;
else
    error('Polyhedral outer approximation not yet implemented');
end

% Input case 1: maxFacets, maxVertices, maxErr
% Input case 2: none


end
