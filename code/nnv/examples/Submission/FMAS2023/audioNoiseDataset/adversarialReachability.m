function [robust,T,PR,T_avg,T_sum] = adversarialReachability(varargin)
%% the function adversarialReachability takes the following inputs
%   nnvnet  : NNV supported neural network model
%   dataset : dataset prepared by the 'createSOCDataset' function
%   percent : the input perturbations as percentage for different noises
%   noiseType   : type of Noise
%   randFeature : the particular feature we want to introduce the noise in

%   and provides the following outputs as a result of approx-star
%   reachability analysis
%   robust  :   local robustness for each of the input audio sequence
%   PR      :   the percentage sample robustness value
%   T_avg   :   the avg time for the reachability calculatiom
%   T_sum   : total time for the reachability calculation

    reachOptions.reachMethod = 'exact-star';
    switch nargin
        case 7
            nnvnet = varargin(1);
            dataset = varargin(2);
            index = varargin(3);
            percent = varargin(4);
            noiseType = varargin(5);
            classIndex = varargin(6);
            randFeature = varargin(7);
        case 6
            nnvnet = varargin(1);
            dataset = varargin(2);
            index = varargin(3);
            percent = varargin(4);
            noiseType = varargin(5);
            classIndex = varargin(6);
        otherwise
            error('Invalid number of inputs, should be 6 or 7');
    end
    percent = percent{1,1};
    nnvnet = nnvnet{1,1};
    randFeature = index{1,1};
    input = dataset{1,1};
    classIndex = classIndex{1,1};
    for j = 1: length(input)
        ip = input(:,:,j)';
        lb = ip;
        ub = ip;
        if strcmp(noiseType{1,1},'SFSI')
            XmuNorm = mean(ip(randFeature,:));
            lb(randFeature,end) = lb(randFeature, end)-XmuNorm*percent;
            ub(randFeature,end) = ub(randFeature, end)+XmuNorm*percent;
        elseif strcmp(noiseType{1,1},'SFAI')
            XmuNorm = mean(ip(randFeature,:));
            lb(randFeature,:) = lb(randFeature,:)-XmuNorm*percent;
            ub(randFeature,:) = ub(randFeature,:)+XmuNorm*percent;
        elseif strcmp(noiseType{1,1},'MFSI')
            XmuNorm = mean(ip(:,:)')';
            lb(:,end) = lb(:, end)-XmuNorm*percent;
            ub(:,end) = ub(:, end)+XmuNorm*percent;
        elseif strcmp(noiseType{1,1},'MFAI')
            XmuNorm = mean(ip(:,:)')';
            lb = lb-XmuNorm*percent;
            ub = ub+XmuNorm*percent;
        end
        IM = ImageStar(lb,ub);
        start_time = tic;
        R(j) = nnvnet.reachSequence(IM, reachOptions);
        T(j) = toc(start_time);
        robust(j) = nnvnet.checkRobust(R(j),classIndex);
    end
    PR = sum(robust==1)/100;
    T_avg = mean(T);
    T_sum = sum(T);   
end