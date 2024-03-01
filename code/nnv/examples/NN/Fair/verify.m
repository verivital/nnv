%% Fairness verification of an NN
% Code based on a few examples
% Data Preprocessing: https://github.com/jonathanlxy/MLproject-UCI-Adult-Income-Classification
% Ideas: https://github.com/LebronX/DeepGemini-public/blob/main/src/German/german_fairness_training.py
% Implementattion/Traihttps://github.com/LebronX/DeepGemini-public/blob/main/src/German/german_fairness_training.pyning

t = tic % track total time for training

%% Read data
Train = csvread('finalset_cleaned_train.csv', 1, 0);
Test  = csvread('finalset_cleaned_test.csv', 1, 0);

%% Neural Network
%{
layers = [
    imageInputLayer(imgSize) % image size = [28 28 1] (Height x Width x Channels)
    flattenLayer;
    
    fullyConnectedLayer(30)
    reluLayer;

    fullyConnectedLayer(30)
    reluLayer

    fullyConnectedLayer(30)
    reluLayer

    fullyConnectedLayer(30)
    reluLayer

    fullyConnectedLayer(numClasses) % 10 = number of classes
    softmaxLayer
    classificationLayer]; 
%}