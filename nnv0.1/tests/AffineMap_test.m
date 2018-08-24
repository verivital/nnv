
% This script is used to test the performance of AffineMap operation in
% Polyhedron class. It turns out that this operation consumes only a very
% small amount of time.
N = 1000;
step = 100;
t = [];
for i=1:step:N
    A = rand(i);
    b = rand(i, 1);
    I = Polyhedron(A, b);
    M = rand(i);
    l = rand(i, 1);
    tic;
    Im = I.affineMap(M);
    Im = Im + l;
    ti = toc;
    t = [t ti];
end

I = 1:step:N;
fig = figure;
plot(I, t, 'b--o');
title('affineMap performance');
xlabel('N - number of constraints');
ylabel('t - time (in seconds)');
saveas(fig, 'affineMapPerformance.pdf');