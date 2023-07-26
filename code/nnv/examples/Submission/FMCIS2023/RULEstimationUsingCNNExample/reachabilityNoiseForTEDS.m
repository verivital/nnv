function [LB,UB,T,PR,POR,T_avg,T_sum] = reachabilityNoiseForTEDS(varargin)
%% the function reachabilityNoiseForTEDS takes the following inputs
%   nnvnet  : NNV supported neural network model
%   dataset : dataset prepared by the 'createTEDSDataset' function
%   percent : the input perturbations as percentage for different noises
%   noiseType   : type of Noise
%   randFeature : the particular feature we want to introduce the noise in

%   and provides the following outputs as a result of approx-star
%   reachability analysis
%   LB      :   the lower bound of the output reachable set
%   UB      :   the upper bound of the output reachable set
%   PR      :   the percentage sample robustness value
%   POR     :   the percentage overlap robustness value
%   T_avg   :   the avg time for the reachability calculatiom
%   T_sum   : total time for the reachability calculation
    reachOptions.reachMethod = 'approx-star';
    switch nargin
        case 6
            nnvnet = varargin(1);
            dataset = varargin(2);
            index = varargin(3);
            percent = varargin(4);
            noiseType = varargin(5);
            randFeature = varargin(6);
        case 5
            nnvnet = varargin(1);
            dataset = varargin(2);
            index = varargin(3);
            percent = varargin(4);
            noiseType = varargin(5);
        otherwise
            error('Invalid number of inputs, should be 3 or 4');
    end
    percent = percent{1,1};
    nnvnet = nnvnet{1,1};
    index = index{1,1};
    input = dataset{1,1}.input{index,1};
    % output = input{1,1}.output;
    allowableUB = dataset{1,1}.allowableUB{index,1};
    allowableLB = dataset{1,1}.allowableLB{index,1};
    for j = 1: length(input)
        ip = input{j,1};
        lb = ip;
        ub = ip;
        if strcmp(noiseType{1,1},'SFSI')
            XmuNorm = mean(ip(randFeature{1,1},:));
            lb(randFeature{1,1},end) = lb(randFeature{1,1}, end)-XmuNorm*percent;
            ub(randFeature{1,1},end) = ub(randFeature{1,1}, end)+XmuNorm*percent;
        elseif strcmp(noiseType{1,1},'SFAI')
            XmuNorm = mean(ip(randFeature{1,1},:));
            lb(randFeature{1,1},:) = lb(randFeature{1,1},:)-XmuNorm*percent;
            ub(randFeature{1,1},:) = ub(randFeature{1,1},:)+XmuNorm*percent;
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
        R = nnvnet.reachSequence(IM, reachOptions);
        T(j) = toc(start_time);
        [L, U]= R.getRanges;
        L_lin = fitLinearLine(L);
        U_lin = fitLinearLine(U);
        LB_actual(j) = L(end);
        UB_actual(j) = U(end);
        LB(j) = L_lin(end);
        UB(j) = U_lin(end);
        LB_allow(j) = allowableLB{j,1}(end)-5;
        UB_allow(j) = allowableUB{j,1}(end)+5;
    end
    % L_lin = fitLinearLine(LB);
    % U_lin = fitLinearLine(UB);
    PR = globalrobustnessSingle(LB,UB,LB_allow,UB_allow);
    POR = globalrobustnessSingle_overlap(LB,UB,LB_allow,UB_allow);
    T_avg = mean(T);
    T_sum = sum(T);   
end
    
function y_lin = fitLinearLine(y)
    x = 1:length(y);
    c = polyfit(x,y,1);
    y_lin = polyval(c,x);
end