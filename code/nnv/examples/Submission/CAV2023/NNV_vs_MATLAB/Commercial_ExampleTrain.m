%% Verification of a fullyconnected network using MATLAB's verification algorithm                       

%% Prepare data and networks

% Load test data
% [XTest,TTest] = digitTest4DArrayData;
% X = XTest(:,:,:,1:10);
% label = TTest(1:10);

% Load networks
load("exampleNet.mat");
load("exampleNetRobust.mat");


% Prepare the network for verification by removing the softmax layer.
net = removeLayers(net,"softmax");
%  When you remove layers from a dlnetwork object, the software returns the network as an uninitialized dlnetwork object. To initialize the network, use the initialize function.
net = initialize(net);
netRobust = removeLayers(netRobust, "softmax");
netRobust = initialize(netRobust);

[XTest,TTest] = digitTest4DArrayData;
% Verification of the whole test set can take a long time. Use a subset of the test data for verification.

numObservations = numel(TTest);
numToVerify = 200;

idx = randi(numObservations,numToVerify,1);
X = XTest(:,:,:,idx);
T = TTest(idx);

% Convert the test images to a dlarray object with the data format "SSCB" (spatial, spatial, channel, batch), which represents image data.
X = dlarray(X,"SSCB");


%% Verify Network Robustness
% To verify the adversarial robustness of a deep learning network, use the verifyNetworkRobustness function. The verifyNetworkRobustness 
% function requires the Deep Learning Toolbox™ Verification Library support package.
% 
% To verify network robustness, the verifyNetworkRobustness function checks that, for all inputs between 
% the specified input bounds, there does not exist an adversarial example. The absence of an adversarial example 
% means that, for all images within the input set defined by the lower and upper input bounds, the predicted class label matches the specified label (usually the true class label).
% 
% For each set of input lower and upper bounds, the function returns one of these values:
% "verified" — The network is robust to adversarial inputs between the specified bounds.
% "violated" — The network is not robust to adversarial inputs between the specified bounds.
% "unproven" — The function cannot prove whether the network is robust to adversarial inputs between the specified bounds.
% 
% Create lower and upper bounds for each of the test images. Verify the network robustness to an input perturbation between –0.02 and 0.02 for each pixel.

perturbation = 0.02;

XLower = X - perturbation;
XUpper = X + perturbation;

% Verify the network robustness for the specified input bounds and true class labels.

result = verifyNetworkRobustness(net,XLower,XUpper,T);
% summary(result)

figure
bar(countcats(result))
xticklabels(categories(result))
ylabel("Number of Observations")


%% Verify Adversarially Trained Network
% Adversarial training is a technique for training a network so that it is robust to adversarial examples [3]. Load a pretrained network 
% that has been trained to be robust to adversarial examples using the methods described in Train Image Classification Network Robust 
% to Adversarial Examples. This network has the same layers as the normal network. The network has been trained to be robust to pixel perturbations in the range [–0.02, 0.02].

load("digitsRobustClassificationMLPNet.mat")
% Prepare the network for verification using the same steps as for the normal network.

netRobust = removeLayers(netRobust,"softmax");
netRobust = initialize(netRobust);


%% Verify the network robustness.

resultRobust = verifyNetworkRobustness(netRobust,XLower,XUpper,T);
% summary(resultRobust)

% Compare the results from the two networks. The robust network has a greater number of observations that correspond to a verified result in comparison to the network without adversarial training.

combineResults = [countcats(result),countcats(resultRobust)];
figure
bar(combineResults)
xticklabels(categories(result))
ylabel("Number of Observations")
legend(["Normal Network","Robust Network"],Location="northwest")


%% Compare Perturbation Values
% Compare the number of verified results as the perturbation value changes. Create lower and upper bounds for each image for a range of perturbation values. To reduce computation time, specify multiple pairs of input bounds in a single call to the verifyNetworkRobustness function.

perturbationRange = 0:0.005:0.05;

XLower = [];
XUpper = [];
TRange = [];

j = 1;
for i = 1:numel(perturbationRange)
    idxRange = j:(j+numToVerify-1);

    perturbationRangeIdx(i,1) = idxRange(1);
    perturbationRangeIdx(i,2) = idxRange(end);

    XLower(:,:,:,idxRange) = X - perturbationRange(i);
    XUpper(:,:,:,idxRange) = X + perturbationRange(i);

    TRange(idxRange) = T;
    j = j + numToVerify;
end

XLower = dlarray(XLower,"SSCB");
XUpper = dlarray(XUpper,"SSCB");

% Verify the robustness of both networks for each pair of input lower and upper bounds.

result = verifyNetworkRobustness(net,XLower,XUpper,TRange);
resultRobust = verifyNetworkRobustness(netRobust,XLower,XUpper,TRange);

% Find the number of verified results for each perturbation value.

numVerified = [];
numVerifiedRobust = [];

for i = 1:numel(perturbationRange)
    range = perturbationRangeIdx(i,:);

    numVerified(i) = sum(result(range(1):range(2)) == "verified");
    numVerifiedRobust(i) = sum(resultRobust(range(1):range(2)) == "verified");
end

% Plot the results. As the perturbation increases, the number of observations returning verified decreases for both networks.
figure
plot(perturbationRange,numVerified,"*-")
hold on
plot(perturbationRange,numVerifiedRobust,"*-")
hold off

legend(["Normal Network","Robust Network"])
xlabel("Perturbation")
ylabel("Number of verified Results")
