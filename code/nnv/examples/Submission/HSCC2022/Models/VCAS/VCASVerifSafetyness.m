function isSafe = VCASVerifSafetyness(init_set)

% Checks if we are in a Safe scenario: verifies whether h remains bounded. 
%
% INPUTS
%
% init_set: state of the physical component of the system in Star form.
%
% OUTPUTS
%
% isSafe: boolean that returns True if we are in a safe state.

isSafe = true;

hlim = 100; % [ft]
epsilon = 10^-3; % to avoid problems with 0
[mh, Mh] = init_set(1).getRange(1); % In theory, init set for VCAS will never be length > 1 because we don't split the angles
[mtau, Mtau] = init_set(1).getRange(4); % In theory, if we don't add unc to tau mtau = Mtau

if  ( (0 <= abs(mtau) && abs(mtau) <= epsilon) || (0 <= abs(Mtau) && abs(Mtau) <= epsilon)) && (abs(mh) < hlim || abs(Mh) < hlim || sign(mh) ~= sign(Mh))
    isSafe = false;
end


end