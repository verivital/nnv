G = [1 0 0; 1 -1 1];
g = [1; 2];

U = HalfSpace(G, g);

x = [2; 1; 0];

bool = U.contains(x);