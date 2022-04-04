function [angSetout] = LimitAngleSet_PiPi(varargin)

% Limits angle set to a -pi to pi interval. Spliting the initial set in two
% set if needed.
% 
% INPUTS
%
% varargin: nargin = 1 --> one star set, angle set in radians.
%           margin = 2 --> (lb,ub) lower and upper bound in radians.
%
% OUTPUTS
%
% angSetout: star set output for angle range between -pi and pi.


if nargin == 1
    inputAngSet = varargin{1};
    if string(class(inputAngSet))~= "Star" || length(inputAngSet) > 1
        error('Wrong input format. Input set must be a single Star set')
    end
    [lb,ub] = inputAngSet.getRanges;
elseif nargin == 2
    lb = varargin{1};
    ub = varargin{2};
end
    % Sets angles in the [-pi,pi] interval
    lb = set_angleRange(lb);
    ub = set_angleRange(ub);
    % Analize if the interval has been divided
    if lb > ub
        % Means that we need to split the set into 2
        angSetout = [Star(-pi,ub) Star(lb,pi)];
    else
        angSetout = Star(lb,ub);
    end
end

