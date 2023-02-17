function Up = PostProcessing(Idx,listofcommands)

% Selects which physical command corresponds to the label 
% (Idx) obtained from the previous step. It's purpose is to correctly store
% all the advisory of all the existing branches at a time step j.
% 
% INPUTS
%
% Idx: index of the min/max value in the weigth vector. Can be a scalar or
%      a vector
% lisofcommands: list of possible commands the controller might take.
%
% OUTPUTS
%
% Up: cell of all of the selected commands in 
%     Star form. Up{i,j} where i is the branch position and j indicates
%     we obtained multiple idx for a branch at the lambda function.


% Set advisory to ownship

Up = {};
for i = 1:length(Idx)
    for j = 1:length(Idx{i})
        Up{i,j} = Star(listofcommands{Idx{i}(j)},listofcommands{Idx{i}(j)});
    end
end

end