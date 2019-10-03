% plot time reduction vs. number of cores 

load reachTime.mat;
load numCores.mat;

fig = figure;
plot(numCores, reachTime, '--*' );
xlabel('N', 'Fontsize', 16);
ylabel('RT', 'Fontsize', 16);
set(gca, 'Fontsize', 16);
saveas(gcf,'time_reduction.pdf');
