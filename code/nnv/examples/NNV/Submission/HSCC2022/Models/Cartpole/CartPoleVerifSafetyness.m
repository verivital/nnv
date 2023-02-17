function isSafe = CartPoleVerifSafetyness(init_set)

% Limits for x and theta found at: https://easychair.org/publications/open/BFKs

% Checks if we are in a Safe scenario: verifies that x, theta and thetadot
% stay bounded . 
%
%
% INPUTS
%
% init_set: state of the physical component of the system in Star form.
%
% OUTPUTS
%
% isSafe: boolean that returns True if we are in a safe state.

isSafe = true;

xlim = 5.0; % [m]
thlim = 24*pi/180; % [rad]
omegalim = 5*pi/180; % [rad] % This one I took it arbitrarly
[mx, Mx] = init_set(1).getRange(1); % In theory, init set for CartPole will never be length > 1 because we don't split the angles
[mth, Mth] = init_set(1).getRange(3);
[momega, Momega] = init_set(1).getRange(4);

if  abs(mth) > thlim || abs(Mth) > thlim 
    isSafe = false;
end


% epsilon = 0.001;
% 
% if ~stop % verifies whether we have violated x's and th's limits 
%     
%     for i = 1:length(step_sets)
%         step_sets(i).getRange(2);
%     end
% end

end