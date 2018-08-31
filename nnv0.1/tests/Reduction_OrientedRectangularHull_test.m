I1 = ExamplePoly.randVrep;   % input set1
I2 = ExamplePoly.randVrep;   % input set2

I = [I1 I2];

V = [0 1; 0 0; 1 0];
I = Polyhedron(V);

P = Reduction.orientedRectangularHull(I);
fig = figure; 
P.plot;
hold on;
I.plot;




