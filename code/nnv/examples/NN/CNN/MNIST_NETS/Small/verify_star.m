load images.mat;
load('Small_ConvNet.mat');
nnvNet = CNN.parse(net, 'Small_ConvNet');

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9


%delta = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
delta = [0.01 0.02 0.03 0.04 0.05];

M = length(delta);
N = 2; % number of images used for robustness tested
del = 5;
inputSetStar = cell(M, 1);
for i=1:M
    inputStar = [];
    for j=1:N
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
        inputStar = [inputStar Z.toImageStar];
    end
    inputSetStar{i} = inputStar;
end

% computing reachable set
numCores = 4;
outputSetStar = cell(M, 1);
verifyTimeStar = zeros(1, M);
for i=1:M
    outputSetStar{i} = nnvNet.reach(inputSetStar{i}, 'approx-star', numCores);
    verifyTimeStar(i) = nnvNet.totalReachTime;
end

%verify results
for i=1:M
    outputSet = outputSetStar{i};
    for j=1:N
        
    end
   
end


