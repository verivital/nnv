function [X, Labels] = sortByLabel(X,Labels)
% SORTBYLABEL Takes an array X which can be anything that corresponds with
% a dlarray dimension ordering. This will sort the array by dlarray
% dimension conventions. 

%   Copyright 2020-2021 The MathWorks, Inc.

Labels = char(Labels); 
numericalLabels = zeros(1, numel(Labels)); 
for i = 1:numel(Labels)
    % SORT the incoming information based on the SCBTU priority.
    switch Labels(i)
        case 'S'
            numericalLabels(i) = 1; 
        case 'C'
            numericalLabels(i) = 2; 
        case 'B'
            numericalLabels(i) = 3; 
        case 'T'
            numericalLabels(i) = 4; 
        case 'U'
            numericalLabels(i) = 5; 
    end
end 

 

% idx should be a stable ordering 
[~, idx] = sort(numericalLabels); 
X = X(idx);
Labels = Labels(idx); 
end
