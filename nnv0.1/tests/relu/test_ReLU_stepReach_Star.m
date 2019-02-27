
V = [1 1 0; 0 1 0; 0 0 1];
C = [1 0; -1 0; 0 1; 0 -1];
d = [1; 1; 1; 1];  % -1 <= a[1] <= 1, -1 <= a[2] <= 2

I1 = Star(V, C, d);
B1 = I1.getBox;
lb = B1.lb;
ub = B1.ub;

index = 2; 
R = ReLU.stepReach_Star(I1, index, lb(index), ub(index));

figure;
I1.plot;
figure;
R(1).plot;
figure;
R(2).plot;