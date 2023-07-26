function dataset = createTEDSDataset(input, output, predictedOutput, seqL, epsilon, type)
%% the function createTEDSDataset prepares the TEDS data for reachability analysis 
% and stores in a struct variable 'dataset'
%   input           : single-featured or multi-featured time series
%   output          : desired output for the dataset
%   preictedOutput  : output predicted by the NN model
%   sqeL            : desired sequence length of each segment of time-series
%   delta           : extent of deviation allowed in output (scale of 0-100)%
%   type            : the noise type

%   dataset : a struct, containing the allowable bounds (i.e., desired
%             output +/- 10 used for RUL estimation dataset), the input, the 
%             desired output and the predicted output
    maxL = max(seqL);
    if type == "cons"
        idxEnd = 151;
    elseif type == "rand"
        dataset.idx = randi(size(input,2) - maxL- 100, 1, 100);
    end
    for j = 1 : length(seqL)
        idx = idxEnd - seqL(j) - 100 + 2: idxEnd;
        for i = 1: 100
            windowInput{i,1} = input(:,idx(i):idx(i)+seqL(j)-1);
            windowOutput{i,1} = output(:,idx(i):idx(i)+seqL(j)-1);
            windowPredOutput{i,1} = predictedOutput(:,idx(i):idx(i)+seqL(j)-1);
        end
        dataset.idx{j,1} = idx;
        dataset.input{j,1} = windowInput;
        dataset.output{j,1} = windowOutput;
        dataset.predoutput{j,1} = windowPredOutput;
    
        dataset.allowableLB{j,1} = cellfun(@(x) max(x-epsilon,0),windowOutput,'uniformoutput',false);
        dataset.allowableUB{j,1} = cellfun(@(x) x+epsilon,windowOutput,'uniformoutput',false);
    end
end