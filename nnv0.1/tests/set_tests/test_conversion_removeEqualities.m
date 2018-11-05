% -1 <= x1 <= 1, x2 = 2
A = [1 0; -1 0]; b = [1; 1]; Ae = [0 1]; be = [2];
P1 = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be);

newP1 = Conversion.removeEqualities(P1);
figure;
P1.plot;
figure;
newP1.plot;