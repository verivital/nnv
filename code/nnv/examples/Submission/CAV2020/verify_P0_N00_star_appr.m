function results = verify_P0_N00_star_appr(N1, N2)

addpath(genpath("../../../engine"));
addpath(genpath("../../../tbxmanager"));
addpath("nnet-mat-files/")
load(['ACASXU_run2a_',num2str(N1),'_',num2str(N2),'_batch_2000.mat']);

% edit property here
P0 = 1;
lb = [55947.69;-3.14;-3.14;1145;0];
ub = [60760;3.14;3.14;1200;60];
unsafe_mat = [-1,0,0,0,0];
unsafe_vec = [-1500];

Layers = [];
n = length(b);
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers = [Layers Li];
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');

Layers = [Layers Ln];
F = FFNNS(Layers);

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

I = Star(lb, ub);

numCores = 4;

[R1, ~] = F.reach(I, 'approx-star'); % approximate reach set using star

normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);

t = tic;
fprintf('\nVerifying over-approximate star reach set...');
R1_norm = R1.affineMap(normalized_mat, normalized_vec); % over-approximate normalized reach set using star
S = R1_norm.intersectHalfSpace(unsafe_mat, unsafe_vec);
check_time = toc(t);

if isempty(S)
    safe = 1;
else
    safe = -1;
end

results.safe = safe;
results.set_number = length(F.outputSet);
results.total_time = check_time + F.totalReachTime;

filename = ['../../../../../../logs/logs_nnv_star_appr/P',num2str(P0),'_N',num2str(N1),num2str(N2),'_star_appr.txt'];
fileID = fopen(filename,'w');
if safe
    fprintf(fileID, 'UNSAT\n');
else
    fprintf(fileID, 'SAT\n');
end
fprintf(fileID,'Number of Output Sets: %d\n', length(F.outputSet));
fprintf(fileID,'Total Time: %f\n', check_time + F.totalReachTime);
fclose(fileID);
