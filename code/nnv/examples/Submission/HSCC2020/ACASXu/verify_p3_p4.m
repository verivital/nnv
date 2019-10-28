clear

P0 = 4; N1 = 5; N2 = 9;
for p0=3:4
    star.safe = []; star.set_num = []; star.time = [];
    star_appr.safe = []; star_appr.set_num = []; star_appr.time = [];
    abs.safe = []; abs.set_num = []; abs.time = [];
    zono.safe = []; zono.set_num = []; zono.time = [];
    
    for n1 = 1:N1
        for n2 = 1:N2
            % %%%%%%%%%%%% exact star
            [results, ~] = verify_P0_N00_star(p0, n1, n2);
            star.safe(end+1,1) = results.safe;
            star.set_num(end+1,1) = results.set_number;
            star.time(end+1,1) = results.total_time;
            
            % %%%%%%%%%%%% approximated star
            [results, ~] = verify_P0_N00_star_appr(p0, n1, n2); 
            star_appr.safe(end+1,1) = results.safe;
            star_appr.set_num(end+1,1) = results.set_number;
            star_appr.time(end+1,1) = results.total_time;  
            
            % %%%%%%%%%%%% abstract domain
            [results, ~] = verify_P0_N00_abs(p0, n1, n2); 
            abs.safe(end+1,1) = results.safe;
            abs.set_num(end+1,1) = results.set_number;
            abs.time(end+1,1) = results.total_time;
            
            % %%%%%%%%%%%%% zonotope
            [results, ~] = verify_P0_N00_zono(p0, n1, n2); 
            zono.safe(end+1,1) = results.safe;
            zono.set_num(end+1,1) = results.set_number;
            zono.time(end+1,1) = results.total_time;
            
        end
    end
    
    switch p0
        case 1
            results_p1.star = star;
            results_p1.star_appr = star_appr;
            results_p1.abs = abs;
            results_p1.zono = zono;
            nnv.p1 = results_p1;
        case 2
            results_p2.star = star;
            results_p2.star_appr = star_appr;
            results_p2.abs = abs;
            results_p2.zono = zono; 
            nnv.p2 = results_p2;
        case 3
            results_p3.star = star;
            results_p3.star_appr = star_appr;
            results_p3.abs = abs;
            results_p3.zono = zono;
            nnv.p3 = results_p3;
        case 4
            results_p4.star = star;
            results_p4.star_appr = star_appr;
            results_p4.abs = abs;
            results_p4.zono = zono;  
            nnv.p4 = results_p4;
    end
    
    
end

save history34.mat

%%  
load others_p34.mat % load results of reluplex, marabou and reluval
property_str = ["p3","p4"];
method_str = {"reluplex";"marabou";"maradnc";"reluval";"zono";"abs";"star";"star_appr"};
for p_temp = property_str
    fprintf("Table of "+p_temp+": \n")
    
    SAT = []; UNSAT=[]; UNK= []; TIMEOUT_1h=[]; TIMEOUT_2h=[]; TIMEOUT_10h=[]; TIME = []; 
    for n = 1:length(method_str)
        method_temp = method_str{n};
        if n<=4
            time_temp = others.(p_temp).(method_temp).time;
            safe_temp = others.(p_temp).(method_temp).safe;
        else
            time_temp = nnv.(p_temp).(method_temp).time;
            safe_temp = nnv.(p_temp).(method_temp).safe;
        end
        
        SAT = [SAT; length(find(safe_temp==0))]; % unsafe
        UNSAT = [UNSAT; length(find(safe_temp==1))]; % safe
        UNK = [UNK; length(find(safe_temp==-1))]; % unknown
        TIMEOUT_1h = [TIMEOUT_1h; length(find(time_temp>=3600*1))];
        TIMEOUT_2h = [TIMEOUT_2h; length(find(time_temp>=3600*2))];
        TIMEOUT_10h = [TIMEOUT_10h; length(find(time_temp>=3600*10))];
        TIME = [TIME; sum(time_temp)]; % total time
        
    end
    T = table(method_str,SAT,UNSAT,UNK,TIMEOUT_1h,TIMEOUT_2h,TIMEOUT_10h, TIME)
    filename = ["Table-of-"+p_temp+".txt"];
    writetable(T,filename)
end


%% helper functions
function [results, R1]= verify_P0_N00_abs(P0, N1, N2)

% addpath(genpath("../../engine"));
% addpath(genpath("../../tbxmanager"));
% addpath("nnet-mat-files/")
load(['ACASXU_run2a_',num2str(N1),'_',num2str(N2),'_batch_2000.mat']);

switch P0
    case 1
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];
        unsafe_mat = [-1 0 0 0 0];
        unsafe_vec = [-1500];
    case 2
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];
        unsafe_mat = [-1 1 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1];
        unsafe_vec = [0; 0; 0; 0];
    case 3
        lb = [1500; -0.06; 3.1; 980; 960];
        ub = [1800; 0.06; 3.14; 1200; 1200];
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
    case 4
        lb = [1500; -0.06; 0; 1000; 700];
        ub = [1800; 0.06; 0; 1200; 800];
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
end

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


[R1, ~] = F.reach(I, 'abs-dom'); % approximate reach set using abstract domain

normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);

t = tic;
R1_norm = R1.affineMap(normalized_mat, normalized_vec); % over-approximate normalized reach set using abstract domain
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
%save(['batch_abs\P',num2str(P0),'_N',num2str(N1),num2str(N2),'_abs.mat'],'results')
end

function [results, R1] = verify_P0_N00_zono(P0, N1, N2)

% addpath(genpath("../../engine"));
% addpath(genpath("../../tbxmanager"));
% addpath("nnet-mat-files/")
load(['ACASXU_run2a_',num2str(N1),'_',num2str(N2),'_batch_2000.mat']);

switch P0
    case 1
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];
        unsafe_mat = [-1 0 0 0 0];
        unsafe_vec = [-1500];
    case 2
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];
        unsafe_mat = [-1 1 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1];
        unsafe_vec = [0; 0; 0; 0];
    case 3
        lb = [1500; -0.06; 3.1; 980; 960];
        ub = [1800; 0.06; 3.14; 1200; 1200];
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
    case 4
        lb = [1500; -0.06; 0; 1000; 700];
        ub = [1800; 0.06; 0; 1200; 800];
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
end

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


[R1, ~] = F.reach(I.getZono, 'approx-zono'); % approximate reach set using zonotope

normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);

t = tic;
R1_norm = R1.affineMap(normalized_mat, normalized_vec); % over-approximate normalized reach set using zonotope
R = R1_norm.toStar;
S = R.intersectHalfSpace(unsafe_mat, unsafe_vec);

check_time = toc(t);
if isempty(S)
    safe = 1;
else
    safe = -1;
end
results.safe = safe;
results.set_number = length(F.outputSet);
results.total_time = check_time + F.totalReachTime;
%save(['batch_zono\P',num2str(P0),'_N',num2str(N1),num2str(N2),'_zono.mat'],'results')

end

function [results, R1] = verify_P0_N00_star_appr(P0, N1, N2)

% addpath(genpath("../../engine"));
% addpath(genpath("../../tbxmanager"));
% addpath("nnet-mat-files/")
load(['ACASXU_run2a_',num2str(N1),'_',num2str(N2),'_batch_2000.mat']);

switch P0
    case 1
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];
        unsafe_mat = [-1 0 0 0 0];
        unsafe_vec = [-1500];
    case 2
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];
        unsafe_mat = [-1 1 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1];
        unsafe_vec = [0; 0; 0; 0];
    case 3
        lb = [1500; -0.06; 3.1; 980; 960];
        ub = [1800; 0.06; 3.14; 1200; 1200];
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
    case 4
        lb = [1500; -0.06; 0; 1000; 700];
        ub = [1800; 0.06; 0; 1200; 800];
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
end

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
%save(['batch_star_appr\P',num2str(P0),'_N',num2str(N1),num2str(N2),'_star_appr.mat'],'results')

end

function [results, R1] = verify_P0_N00_star(P0, N1, N2)

% addpath(genpath("../../engine"));
% addpath(genpath("../../tbxmanager"));
% addpath("nnet-mat-files/")
load(['ACASXU_run2a_',num2str(N1),'_',num2str(N2),'_batch_2000.mat']);

switch P0
    case 1
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];
        unsafe_mat = [-1 0 0 0 0];
        unsafe_vec = [-1500];
    case 2
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];
        unsafe_mat = [-1 1 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1];
        unsafe_vec = [0; 0; 0; 0];
    case 3
        lb = [1500; -0.06; 3.1; 980; 960];
        ub = [1800; 0.06; 3.14; 1200; 1200];
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
    case 4
        lb = [1500; -0.06; 0; 1000; 700];
        ub = [1800; 0.06; 0; 1200; 800];
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
end

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

c = parcluster('local');
numCores = c.NumWorkers;

[R1, ~] = F.reach(I, 'exact-star', numCores); % exact reach set using star

normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);

t = tic;
n = length(R1);
R1_norm = [];
parfor i=1:n
    R1_norm = [R1_norm  R1(i).affineMap(normalized_mat, normalized_vec)]; % exact normalized reach set
end
check_time = toc(t);

t = tic;
fprintf('\nVerifying exact star reach set...');
unsafe = 0;
n = length(R1);
parfor i=1:n
    S = R1_norm(i).intersectHalfSpace(unsafe_mat, unsafe_vec);
    if ~isempty(S)
        unsafe = unsafe + 1;
    end
end

check_time = check_time + toc(t);

if unsafe>=1
    safe = 0;
else
    safe = 1;
end

results.safe = safe;
results.set_number = length(F.outputSet);
results.total_time = check_time + F.totalReachTime;
%save(['batch_star\P',num2str(P0),'_N',num2str(N1),num2str(N2),'_star.mat'],'results')

end

