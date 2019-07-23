function obj = Zonotope(c, G)
% Constructs a zonotope based on the center point and generators
%
% Z = Zonotope(C, G) constructs a zonotope with the center
% point C and the generators G. The center C must be a n-by-1
% vector and the generators must be n-by-m matrices where "m"
% is the number of generators.

narginchk(2, 2);
assert(size(c, 2)==1, 'The center point must be a column vector.');
assert(size(G, 1)==numel(c), 'The generators must be a set of %dx1 vectors.', numel(c));

% compute vertices by adding the generators to the center point
V = c;
for i = 1:size(G, 2)
    V = [bsxfun(@plus, V, G(:, i)), bsxfun(@plus, V, -G(:, i))];
end
obj = Polyhedron(V'); % in Polyhedron, vertices are stored row-wise

end
