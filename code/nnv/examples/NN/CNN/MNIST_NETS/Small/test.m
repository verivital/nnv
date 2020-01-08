load images.mat;
load('Small_ConvNet.mat');

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9


nnvNet = CNN.parse(net, 'Small_ConvNet');
IM = IM_data(:,:,1);
IM = reshape(IM, [28, 28, 1]);
N = 784; 
IM1 = reshape(IM, [N 1]);

% Brightening attack
del1 = 15;
lb = IM1;
ub = IM1;
for i=1:N
    if  IM1(i) >= 255 - del1
        lb(i) = 0.3*IM1(i);
        ub(i) = 0.5*IM1(i);
    end
end

lb = reshape(lb, [28 28 1]);
ub = reshape(ub, [28 28 1]);

inputSet = ImageZono(lb, ub);

figure;
imshow(uint8(IM));
figure;
imshow(uint8(lb));
figure;
imshow(uint8(ub));

y1 = nnvNet.evaluate(lb);

%outputSet = nnvNet.reach(inputSet, 'approx-zono');
