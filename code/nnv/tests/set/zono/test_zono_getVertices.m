c = [1; 1];
V = [2 1 1; -1 1 0];
Z = Zono(c, V);
V = Z.getVertices;

P = Polyhedron('V', V');

figure;
Z.plot;
figure;
P.plot;