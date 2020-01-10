load images.mat;
load('Small_ConvNet.mat');
nnvNet = CNN.parse(net, 'Small_ConvNet');

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9


%delta = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
delta = [0.005 0.01 0.015 0.02 0.025];
%delta = 0.01;
M = length(delta);
N = 100; % number of test images used for testing robustness
del = 5;
inputSetStar = cell(M, 1);
for i=1:M
    inputStar = [];
    count = 0;
    for j=1:2000
        IM = IM_data(:,:,j);
        lb = IM;
        ub = IM;
        for k=1:784
            if  IM(k) >= 255 - del
                lb(k) = 0;
                ub(k) = delta(i)*IM(k);
            end
        end
        lb = reshape(lb, [28 28 1]);
        ub = reshape(ub, [28 28 1]);
        Z = ImageZono(lb, ub);
        S = Z.toImageStar;
        if ~isempty(S.C)
            inputStar = [inputStar S];
            count = count + 1;
        end
        if count == N
            break;
        end
    end
    inputSetStar{i} = inputStar;
end

% evaluate robustness
numCores = 4;
verifyTimeAbsDom = zeros(1, M);
correct_labels = IM_labels(1:N, 1);
r_abs_dom = zeros(1, M); % robustness percentage on an array of test N input sets
ids_abs_dom = cell(1,M); % classified labels
for i=1:M
    t = tic;
    [r_abs_dom(i),ids_abs_dom{i}] = nnvNet.evaluateRobustness(inputSetStar{i}, correct_labels, 'abs-dom', numCores);
    verifyTimeAbsDom(i) = toc(t);
end

save result_abs_dom.mat verifyTimeAbsDom r_abs_dom ids_abs_dom


