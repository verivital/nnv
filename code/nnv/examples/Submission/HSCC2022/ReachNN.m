function y = ReachNN(Unn,listNN,method,combos,agent)

% Choose NN to execute and computes the NN algebraic operations.
% 
% INPUTS
%
% UNN: array with the normalized state variables of every possible branch
%      in star form.
% listNN: list with the info of all NNs used.
% method:reach method used by the neural network.
% combos: vector containing prev advisory, output set, current advisory,
%         state set.
% agent: id of the agent
%
% OUTPUTS
%
% y: cell array with the weighted vector of all possible branches.

y = {};

for i = 1:size(Unn,2)
    
    % Indexes of the branches where Output set is the same 
    idxs = find(cell2mat(combos(:,2)) == i);
    idx = idxs(1);
    % Choice of NN
    prevadv = combos{idx,1}(agent);
    nn_name = sprintf('n%d',prevadv); % If we don't name n1,n2 etc the diferent NN it is not going to work beacuse of syntaxis
    % Computation of such NN
    ya = listNN.(nn_name).reach(Unn(i),method);
    y{i} = ya;
    
end

end
