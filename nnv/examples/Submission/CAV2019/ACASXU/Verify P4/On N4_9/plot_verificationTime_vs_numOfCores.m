load safety_checking_time.mat;
load reachTime.mat;
load numCores.mat;

n = length(reachTime);

verificationTime = zeros(n, 1);
ReluplexTime = zeros(n,1);
for i=1:n
    verificationTime(i) = reachTime(i) + safetycheckingTime(i);
    ReluplexTime(i) = 489;
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