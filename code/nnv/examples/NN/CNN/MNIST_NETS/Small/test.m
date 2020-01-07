load images.mat;
load('Small_ConvNet.mat');

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9


nnvNet = CNN.parse(net, 'Small_ConvNet');
IM = IM_data(:,:,1);
IM = reshape(IM, [28, 28, 1]);
imshow(IM);

N = 784; 
IM1 = reshape(IM, [N 1]);

% Brightening attack
del = 50;
lb = zeros(N, 1);
ub = zeros(N, 1);
for i=1:n
    if IM1(i) >= 255 - del
        lb(i) = IM1(i);
        ub(i) = 255;
    else
        lb(i) = IM1(i);
        ub(i) = IM1(i);
    end
end

lb = reshape(lb, [28 28 1]);
ub = reshape(ub, [28 28 1]);

inputSet = ImageZono(lb, ub);

figure;
imshow(IM);
figure;
imshow(lb);
figure;
imshow(ub);

%outputSet = nnvNet.reach(inputSet, 'approx-zono');
