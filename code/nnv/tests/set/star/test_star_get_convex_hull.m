
I1 = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I1 = Star(V', I1.A, I1.b); % star set 1

I2 = ExamplePoly.randVrep;  
V = [1 1; 1 0; 0 1];
I2 = Star(V', I2.A, I2.b); % star set 2

S = Star.get_convex_hull([I1, I2]);

figure;
S.plot; 
hold on;
I1.plot;
hold on;
I2.plot;

