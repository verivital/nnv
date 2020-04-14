P1 = Polyhedron('A', R1.C, 'b', R1.d);
W = R1.V(:, 2:R1.nVar + 1);
P2 = P1.affineMap(W);