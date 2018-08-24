% This script is used to test the performance of outerApprox operation in
% mpt toolbox. It turns out that this operation 

N = 200;
step = 10;
t = [];
for i=1:step:N
    A = rand(i);
    b = rand(i, 1);
    I = Polyhedron(A, b);
    tic;
    outerApprox(I);
    ti = toc;
    t = [t ti];
end

I = 1:step:N;
fig = figure;
plot(I, t, 'b--o');
title('outerApprox performance');
xlabel('N - number of constraints');
ylabel('t - time (in seconds)');
saveas(fig, 'outerApproxPerformance.pdf');