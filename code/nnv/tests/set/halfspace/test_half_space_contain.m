% HalfSpace class defining Gx <= g

% Halfspace define by the region 
%  -x1 <= 5
%   x2 <= 0
G = [-1 0; 0 1];
g = [5; 0];

U = HalfSpace(G, g);

x = [2; -1];

bool = U.contains(x);
assert(bool == 1)

x1 = [-2; 1];

bool1 = U.contains(x1);
assert(bool1 == 0);