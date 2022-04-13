function init_set = InitializeState(Param, init_state, unc)

% Creates the initial set in Star form by adding uncertainties to a given 
% initial state.
% 
% INPUTS
%
% Param: struct containing the information needed to execute the simulation.
% init_state: a real-valued vector representing an initial state of the
% physical component of the system.
% unc: a real-valued vector representing the uncertainties on the initial
% state of the physical component of the system.
%
% OUTPUTS
%
% init_set: initial set in Star form. 

geom = Param.geometry;
lb = init_state';
ub = init_state';
dim = length(init_state);

for k=1:dim
    lb(k) = lb(k) - unc(k);
    ub(k) = ub(k) + unc(k);
end

% Converts intervals in Star object
init_set = geom(lb,ub);

end