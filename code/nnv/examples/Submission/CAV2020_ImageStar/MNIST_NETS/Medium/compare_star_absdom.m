load images.mat;
load('Medium_ConvNet.mat');
nnvNet = CNN.parse(net, 'Medium_ConvNet');

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9



delta = [0.005 0.01 0.015];
d = [250 245 240]; % threshold for brightening attack

% for testing
%delta = 0.01; 
%d = 240;

M = length(delta);
P = length(d);
N = 100; % number of test images used for testing robustness

%% construct input sets with different values of d and delta

inputSetStar = cell(P, M);
correct_labels = cell(P, M);

for i=1:P
    for j=1:M
        
        inputStar = [];
        count = 0;
        labels = zeros(1, N);
        for k=1:2000
            IM = IM_data(:,:, k);
            lb = IM;
            ub = IM;
            for l=1:784
                if IM(l) >= d(i)
                    lb(l) = 0;
                    ub(l) = delta(j)*IM(l);
                end
            end
            
            lb = reshape(lb, [28 28 1]);
            ub = reshape(ub, [28 28 1]);
            Z = ImageZono(lb, ub);
            S = Z.toImageStar;
            if ~isempty(S.C)
                inputStar = [inputStar S];
                count = count + 1;
                labels(count) = IM_labels(k);
            end            
            if count == N
                break;
            end
            
        end        
        inputSetStar{i, j} = inputStar;
        correct_labels{i, j} = labels;
    end
end



%% evaluate robustness

VT_star = zeros(P, M); % verification time of the approx-star method
VT_absdom = zeros(P, M); % verification time of the DeepPoly abstract domain method

r_star = zeros(P, M); % robustness percentage on an array of N tested input sets obtained by the approx-star method
r_absdom = zeros(P, M); % robustness percentage on an array of N tested input sets obtained by the DeepPoly abstract domain method

c = parcluster('local');
numCores = c.NumWorkers; % specify number of cores used for verification


for i=1:P
    for j=1:M
             
        t = tic;
        r_star(i, j) = nnvNet.evaluateRobustness(inputSetStar{i, j}, correct_labels{i, j}, 'approx-star', numCores);
        VT_star(i, j) = toc(t);
                
        t = tic;
        r_absdom(i, j) = nnvNet.evaluateRobustness(inputSetStar{i, j}, correct_labels{i, j}, 'abs-dom', numCores);
        VT_absdom(i, j) = toc(t);
        
        
    end
end

save Medium_ConvNet_Results.mat r_star VT_star r_absdom VT_absdom;


%% print the results

fprintf('\n========================================================================================');
fprintf('\n          ROBUSTNESS VERIFICATION RESULTS (IN PERCENT) OF MEDIUM_CONVNET                 ');
fprintf('\n========================================================================================\n\n');


for j=1:M
    fprintf("             delta = %.5f", delta(j));
end

fprintf("\n");

for j=1:M
    fprintf("         Polytope   ImageStar");
end

fprintf("\n");

for i=1:P
    fprintf("d = %d", d(i));
    for j=1:M
        fprintf("    %.2f        %.2f      ", 100*r_absdom(i, j), 100*r_star(i, j));
    end
    fprintf("\n");
end

fprintf('\n========================================================================================');
fprintf('\n                VERIFICATION TIMES (IN SECONDS) OF MEDIUM_CONVNET                        ');
fprintf('\n========================================================================================\n\n');


for j=1:M
    fprintf("             delta = %.5f", delta(j));
end

fprintf("\n");

for j=1:M
    fprintf("         Polytope   ImageStar");
end

fprintf("\n");

for i=1:P
    fprintf("d = %d", d(i));
    for j=1:M
        fprintf("    %.2f       %.2f        ", VT_absdom(i, j), VT_star(i, j));
    end
    fprintf("\n");
end


%% Print to file
fid = fopen('Medium_ConvNet_Results.txt', 'wt');
fprintf(fid,'\n========================================================================================');
fprintf(fid,'\n          ROBUSTNESS VERIFICATION RESULTS (IN PERCENT) OF MEDIUM_CONVNET                 ');
fprintf(fid,'\n========================================================================================\n\n');


for j=1:M
    fprintf(fid,"             delta = %.5f", delta(j));
end

fprintf(fid,"\n");

for j=1:M
    fprintf(fid,"         Polytope   ImageStar");
end

fprintf(fid,"\n");

for i=1:P
    fprintf(fid,"d = %d", d(i));
    for j=1:M
        fprintf(fid,"    %.2f        %.2f      ", 100*r_absdom(i, j), 100*r_star(i, j));
    end
    fprintf(fid,"\n");
end

fprintf(fid,'\n========================================================================================');
fprintf(fid,'\n                VERIFICATION TIMES (IN SECONDS) OF Medium_CONVNET                        ');
fprintf(fid,'\n========================================================================================\n\n');


for j=1:M
    fprintf(fid,"             delta = %.5f", delta(j));
end

fprintf(fid,"\n");

for j=1:M
    fprintf(fid,"         Polytope   ImageStar");
end

fprintf(fid,"\n");

for i=1:P
    fprintf(fid,"d = %d", d(i));
    for j=1:M
        fprintf(fid,"    %.2f       %.2f        ", VT_absdom(i, j), VT_star(i, j));
    end
    fprintf(fid,"\n");
end

%% Print latex table
fid = fopen('Medium_ConvNet_Results.tex', 'wt');
fprintf(fid,'\nRobustness results\n');
for i=1:P
    fprintf(fid, '$d = %d$', d(i));
    for j=1:M
        fprintf(fid, '& $%.2f$ & $%.2f$', 100*r_absdom(i, j), 100*r_star(i,j));
    end
    fprintf(fid,"\\\\");
    fprintf(fid,"\n");
end

fprintf(fid,'\nVerification Times\n');
for i=1:P
    fprintf(fid, '$d = %d$', d(i));
    for j=1:M
        fprintf(fid, '& $%.2f$ & $%.2f$', VT_absdom(i, j), VT_star(i,j));
    end
    fprintf(fid,"\\\\");
    fprintf(fid,"\n");
end
