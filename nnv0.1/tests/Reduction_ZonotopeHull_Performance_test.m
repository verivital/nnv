
I = ExamplePoly.randVrep;   % input set
N = 20;
step = 5;
t = [];
for i=1:step:N
    A = rand(i, 2);
    b = rand(i, 1);
    P = I.affineMap(A);
    P = P + b;
    P.minVRep();

    tic;
    Reduction.zonotopeHull_Random(P.V', 10);
    ti = toc;
    t = [t ti];
    
end

I = 1:step:N;
fig = figure;
plot(I, t, 'b--o');
title('zonotopeHull_Random performance');
xlabel('N - number of constraints');
ylabel('t - time (in seconds)');
saveas(fig, 'zonotopeHullPerformance.pdf');