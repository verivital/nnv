load safety_checking_time.mat;
load reachTime.mat;
load numCores.mat;

n = length(reachTime_star);

verificationTime = zeros(n, 1);
ReluplexTime = zeros(n,1);
for i=1:n
    verificationTime(i) = reachTime_star(i) + safetycheckingTime(i);
    ReluplexTime(i) = 653;
end

fig = figure;
plot(numCores, verificationTime, '-*');
hold on;
plot(numCores, ReluplexTime, '-x');
xlabel('N', 'Fontsize', 16);
ylabel('VT (sec)', 'Fontsize', 16);
legend('Star', 'Reluplex');
set(gca, 'Fontsize', 16);
saveas(gcf, 'verificationTime_vs_numCores.pdf');