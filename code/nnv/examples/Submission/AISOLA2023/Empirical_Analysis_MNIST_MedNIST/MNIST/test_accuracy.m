%% We are going to analyze the test accuracy across classes here

%% 1) Load data 

% Load data
dataPath = "../../../../../../../data/MNIST/";
filenameImagesTest = 't10k-images-idx3-ubyte.gz';
filenameLabelsTest = 't10k-labels-idx1-ubyte.gz';

% use test data for verification
XData = processImagesMNIST(dataPath + filenameImagesTest);
XData_dl = dlarray(XData, "SSCB"); % for dlnetworks
YData = processLabelsMNIST(dataPath + filenameLabelsTest);
YData = double(YData)-1; % MATLAB network outputs goes from 0 to 9, but data goes from 1 to 10;


%% 2) Evaluate all models on the complete dataset

% Iterate trhough every folder and subfolder and analyze all the models
path = pwd;
folders = dir(path);
% Skip the first two that appear in every folder and subfolder as those
% correspond to (".", and "..")

N = 10000;
results = zeros(45,N); % to save the predicted results
model_count = 1;

accuracy = zeros(45,1);
models(45) = string;

% Go into every folder of and analyze each model
for r = 4:6 % iterate through regularizers (3)
    sub_path = [path, filesep, folders(r).name, filesep];
    inits_path = dir(sub_path);
    for i = 3:length(inits_path) % go through all initializations (3 x 3)
        if inits_path(i).isdir
            temp_path = [sub_path, inits_path(i).name, filesep, 'models', filesep];
            models_path = dir([temp_path, '*.mat']);
            for m = 1:length(models_path) % go through all models ( 5 x 3 x 3 )
                netpath = [temp_path, models_path(m).name];
                net_info = load(netpath);
                curr_net = net_info.net;
                if isa(curr_net, "dlnetwork")
                    yPred = predict(curr_net, XData_dl);
                    [~,predictedLabels] = max(extractdata(yPred));
                    predictedLabels = predictedLabels'-1;
                else
                    predictedLabels = classify(curr_net, XData);
                    predictedLabels = cellfun(@(x) str2double(x), cellstr(predictedLabels));
                end
                acc_labels = predictedLabels == YData;
                accuracy(model_count) = sum(acc_labels)/N;
                results(model_count,:) = acc_labels;
                models(model_count) = models_path(m).name;
                model_count = model_count + 1;
            end
        end
    end
end
       
   
% data indexes (balanced data, but not same number of images for all)
Zero  = find(YData == 0);
One   = find(YData == 1);
Two   = find(YData == 2);
Three = find(YData == 3);
Four  = find(YData == 4);
Five  = find(YData == 5);
Six   = find(YData == 6);
Seven = find(YData == 7);
Eight = find(YData == 8);
Nine  = find(YData == 9);

acc_classes = zeros(45,10);
for i = 1:45
    acc_classes(i,1)  = sum(results(i,Zero)) / length(Zero);
    acc_classes(i,2)  = sum(results(i,One))  / length(One);
    acc_classes(i,3)  = sum(results(i,Two))  / length(Two);
    acc_classes(i,4)  = sum(results(i,Three))/ length(Three);
    acc_classes(i,5)  = sum(results(i,Four)) / length(Four);
    acc_classes(i,6)  = sum(results(i,Five)) / length(Five);
    acc_classes(i,7)  = sum(results(i,Six))  / length(Six);
    acc_classes(i,8)  = sum(results(i,Seven))/ length(Seven);
    acc_classes(i,9)  = sum(results(i,Eight))/ length(Eight);
    acc_classes(i,10) = sum(results(i,Nine)) / length(Nine);
end



%% Finally, save images that all models correctly predict to analyze its robustness

% we are going to save the first 50 indexes per class that are correctly
% classified by all models
numClasses = 10;
nClassImgs = 30;

xVerIdxs = zeros(numClasses*nClassImgs,1);
class_idxs = {Zero, One, Two, Three, Four, Five, Six, Seven, Eight, Nine};

for cl=1:10
    xVerIdxs((cl-1)*nClassImgs+1:cl*nClassImgs) = findAccData(results, class_idxs{cl}, nClassImgs);
end


% Save results and verification image indexes
save("acc_results.mat", "xVerIdxs", "results", "acc_classes", "accuracy", "models");


%% Helper Functions
function good_idxs = findAccData(res_models,idxs, nVerif)
    k = 1;
    good_idxs = zeros(nVerif,1);
    for ix = idxs'
        if all(res_models(:,ix) == 1)
            good_idxs(k) = ix;
            k = k+1;
        end
        if k > nVerif
            break
        end
    end
end
