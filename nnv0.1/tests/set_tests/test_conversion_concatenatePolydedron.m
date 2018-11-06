% -1 <= x1 <= 1, x2 = 2

P1 = Polyhedron('lb', [-1; 0], 'ub', [1; 2]);

P2 = Polyhedron('lb', 0, 'ub', 1);

P = Conversion.concatenatePolyhedron([P1 P2]);

figure; 
P.plot;


