
I1 = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I1 = Star(V', I1.A, I1.b); % star set 1

I2 = ExamplePoly.randVrep;  
V = [1 1; 1 0; 0 1];
I2 = Star(V', I2.A, I2.b); % star set 2


I3 = ExamplePoly.randVrep;  
V = [1 2; 1 0; 1 1];
I3 = Star(V', I3.A, I3.b); % star set 3


I4 = ExamplePoly.randVrep;  
V = [1 1; 0 4; 3 1];
I4 = Star(V', I4.A, I4.b); % star set 4



I = [I1 I2 I3 I4];

merge1 = Star.merge_stars(I, 2, 'single');


figure;
merge1(1).plot; 

figure; 
merge1(2).plot; 

figure;
I1.plot;
figure;
I2.plot;
figure;
I3.plot;
figure;
I4.plot;


