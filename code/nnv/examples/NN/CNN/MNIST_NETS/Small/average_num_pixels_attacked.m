load images.mat;
% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9



delta = [0.005 0.01 0.015];
d = [250 245 240]; % threshold for brightening attack

% for testing
%delta = 0.005; 
%d = 250;

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

%% Compute average number of pixels attacked
num_pixels_attacked = zeros(P, M);

for i=1:P
    for j=1:M
        count = 0;
        inputSet = inputSetStar{i, j};
        for k=1:N
            count = count + inputSet(k).numPred;
        end
        num_pixels_attacked(i, j) = ceil(count/N);
    end
end

%% print the results

fprintf('\n========================================================================================');
fprintf('\n                             AVERAGE NUMBER OF PIXELS ATTACKED                          ');
fprintf('\n========================================================================================\n\n');


for j=1:M
    fprintf("             delta = %.5f", delta(j));
end

fprintf("\n");

for i=1:P
    fprintf("d = %d", d(i));
    for j=1:M 
        fprintf("             %d           ", num_pixels_attacked(i, j));
    end
    fprintf("\n");
end

