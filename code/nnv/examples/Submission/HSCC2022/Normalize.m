function UNN = Normalize(Ro,mean,range,cmbs,agent)

% Normalizes the input of the NN.
% 
% INPUTS
%
% Ro: array with the state variables of every possible branch in star form.
% mean: array of dimension = number of NN
% range: array of dimension = number of NN
%
% OUTPUTS
%
% UNN: array with the normalized state variables of every possible branch
%      in star form.

% Normalize inputs
dim = size(mean,2);
UNN = [];
for i=1:length(Ro)
    % Normalize rho, theta and psi
    pos = cmbs{i,1}(agent); % previous adv (determines what NN to use)
    Un = Ro(i).affineMap(eye(dim),-1*mean(pos,1:dim)');
    Un = Un.affineMap(diag(1./range(pos,1:dim)),[]);
    UNN = [UNN Un];
end
end