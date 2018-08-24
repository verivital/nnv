% This script is used to test the performance of ConvexHull operation in
% mpt toolbox. It turns out that this operation is very
% time-consumming. One can test and verify this by changing the values of N
% and step in the following script.

N = 10;
step = 1;
t = [];
for i=1:step:N
    A1 = rand(i);
    b1 = rand(i, 1);
    P(1) = Polyhedron(A1, b1);
    A2 = rand(i);
    b2 = rand(i, 1);
    P(2) = Polyhedron(A2, b2);
    U = PolyUnion('Set', P);
    tic;
    U.convexHull;
    ti = toc;
    t = [t ti];
end

I = 1:step:N;
fig = figure;
plot(I, t, 'b--o');
title('convexHull performance');
xlabel('N - number of constraints');
ylabel('t - time (in seconds)');
saveas(fig, 'convexHullPerformance.pdf');