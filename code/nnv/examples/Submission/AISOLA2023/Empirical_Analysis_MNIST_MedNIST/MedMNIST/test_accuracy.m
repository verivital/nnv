%% We are going to analyze the test accuracy across classes here

%% 1) Load data 

dataFolder = "../../../../../../../data/MedNIST";

% Go through every folder (label) and load all images
categs = dir(dataFolder);
% Initialize vars
N = 58954; % total number of images
XData = zeros(64, 64, 1, N); % height = width = 64, 10k imgs per class (greyscale)
YData = zeros(N, 1);        % except BreatMRI, that has only 8954
% Load images
count = 1;
for i = 3:length(categs)-1
    label = dir(dataFolder + "/"+ string(categs(i).name));
    for k = 3:length(label)
        XData(:, :, :, count) = imread([label(k).folder '/' label(k).name]);
        YData(count) = i-2;
        count = count + 1;
    end
end
% YData = categorical(YData);
XData_dl = dlarray(XData, "SSCB");

%% 2) Evaluate all models on the complete dataset

% Iterate trhough every folder and subfolder and analyze all the models
path = pwd;
folders = dir(path);
% Skip the first two that appear in every folder and subfolder as those
% correspond to (".", and "..")

results = zeros(45,58954); % to save the predicted results
model_count = 1;

accuracy = zeros(45,1);
models(45) = string;

% Go into every folder of and analyze each model
for r = 3:5 % iterate through regularizers (3)
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
                    predictedLabels = predictedLabels';
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

% 3) Get accuracy statistics per model per class
% 6 classes: [Abdominal CT, Breast MRI,  Chest CT,   cxr ,  Hand,  Head CT]
%   labels : [      1          2             3        4       5      6   ]           
   
% data indexes
abdomen = 1:10000;
breast = 10001:18954; 
chest = 18955:28954;
cxr = 28955:38954;
hand = 38955:48954;
head = 48955:58954;

acc_classes = zeros(45,6);
for i = 1:45
    acc_classes(i,1) = sum(results(i,abdomen));
    acc_classes(i,2) = sum(results(i,breast));
    acc_classes(i,3) = sum(results(i,chest));
    acc_classes(i,4) = sum(results(i,cxr));
    acc_classes(i,5) = sum(results(i,hand));
    acc_classes(i,6) = sum(results(i,head));
end



%% Finally, save images that all models correctly predict to analyze its robustness

% we are going to save the first 50 indexes per class that are correctly
% classified by all models
numClasses = 6;
nClassImgs = 50;
k = 1;
xVerIdxs = zeros(numClasses*nClassImgs,1);

% 1) Abdomen
xVerIdxs(1:nClassImgs) = findAccData(results, abdomen, nClassImgs);

% 2) Breast
xVerIdxs(nClassImgs + 1:nClassImgs*2) = findAccData(results, breast, nClassImgs);

% 3) chest
xVerIdxs(nClassImgs*2+1:nClassImgs*3) = findAccData(results, chest, nClassImgs);

% 4) cxr
xVerIdxs(nClassImgs*3+1:nClassImgs*4) = findAccData(results, cxr, nClassImgs);

% 5) Hand
xVerIdxs(nClassImgs*4+1:nClassImgs*5) = findAccData(results, hand, nClassImgs);

% 6) Head
xVerIdxs(nClassImgs*5+1:nClassImgs*6) = findAccData(results, head, nClassImgs);


% Save results and verification image indexes
save("acc_results.mat", "xVerIdxs", "results", "acc_classes", "accuracy", "models");


%% Helper Functions
function good_idxs = findAccData(res_models,idxs, nVerif)
    k = 1;
    ix = idxs(1);
    good_idxs = zeros(nVerif,1);
    while k <= nVerif && ix <= idxs(end)
        if all(res_models(:,ix) == 1)
            good_idxs(k) = ix;
            k = k+1;
        end
        ix = ix+1;
    end
end
