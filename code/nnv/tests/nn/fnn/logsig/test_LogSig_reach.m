
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

B = I.getBox; 
I_Zono = B.toZono;



S1 = LogSig.reach(I, 'approx-star');
S2 = LogSig.reach(I, 'abs-dom');
S3 = LogSig.reach(I_Zono, 'approx-zono');

figure;
I.plot; % input set

figure;
Zono.plots(S3);
hold on;
Star.plots(S2); 
hold on;
Star.plots(S1);
