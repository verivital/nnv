I1 = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I1 = Star(V', I1.A, I1.b); % star set 1

I2 = ExamplePoly.randVrep;  
V = [1 1; 1 0; 0 1];
I2 = Star(V', I2.A, I2.b); % star set 2


I3 = Star.get_hypercube_hull([I1 I2]);

figure;
Box.plot(I3); 
hold on;
Star.plot(I2);
hold on;
Star.plot(I1);