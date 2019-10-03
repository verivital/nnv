clc;
clear;
numSplits = [1 2 3 4 5];
n = length(numSplits);
memUsage = zeros(n,1);
for i=1:n
    memUsage(i) = 224*224*64*8*(numSplits(i))/1024/1024;
end


figure;
plot(numSplits, memUsage, '--x');
ax = gca;
ax.FontSize = 13; 
xlabel('Number of Splits','FontSize',13);
ylabel('Memory Usage (GB)', 'FontSize', 13);



