load m2nist_net_results_attackedpixels.mat;
N1 = Nets(1);
N2 = Nets(2);

N1_relu_reachTime = N1.reachTime(4) + N1.reachTime(8) + N1.reachTime(12);
N1_others_reachTime = sum(N1.reachTime) - N1_relu_reachTime; 

relu_id = [3 5 8 10 13 15 17 19];
N2_relu_reachTime = 0;
for i=1:length(relu_id)
    N2_relu_reachTime = N2_relu_reachTime + N2.reachTime(relu_id(i));
end

N2_other_reachTime = sum(N2.reachTime) - N2_relu_reachTime;

fig = figure; 
Y = [N1_relu_reachTime N1_others_reachTime; N2_relu_reachTime N2_other_reachTime];

subplot(2, 2, 1);
b1 = bar(1:length(N1.reachTime), diag(N1.reachTime), 'stacked', 'FaceColor', 'c');
b1(4).FaceColor = 'r';
b1(8).FaceColor = 'r';
b1(12).FaceColor = 'r';

% xtips1 = [4 8 12];
% ytips1 = N1.reachTime(xtips1);
% labels1 = 'R';
% text(xtips1, ytips1, labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
ylim([0 max(N1.reachTime)+10]);
xlabel('Layer ID (N_4)');
ylabel('Reach-Time (sec)');

subplot(2, 2, [2 4]);
b3 = bar(Y);
b3(1).FaceColor = 'r';
b3(2).FaceColor = 'c';
set(gca, 'XTickLabel', {'N_4', 'N_5'});
xlabel('Network');
ylabel('Total Reach-Time (sec)');


subplot(2, 2, 3);
b2 = bar(1:length(N2.reachTime), diag(N2.reachTime), 'stacked', 'FaceColor', 'c');
for i=1:length(relu_id)
    b2(relu_id(i)).FaceColor = 'r';
end
% xtips2 = relu_id;
% ytips2 = N2.reachTime(xtips2);
% labels2 = 'R';
% text(xtips2, ytips2, labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
ylim([0 max(N2.reachTime)+20]);
xlabel('Layer ID (N_5)');
ylabel('Reach-Time (sec)');