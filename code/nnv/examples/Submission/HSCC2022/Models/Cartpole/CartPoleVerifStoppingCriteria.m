function stop = CartPoleVerifStoppingCriteria(isSafe, init_set)

% Stops the simulation under certain conditions.
%
% INPUTS
%
% isSafe: boolean that returns True if we are in a safe state. 
% init_set: state of the physical component of the system in Star form.
%
% OUTPUTS
%
% stop: boolean than returns True if the simulations should be stopped.

stop = ~isSafe; 

end