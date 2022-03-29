function u = ACASXuCommands

% Returns the possible commands that the controller can execute.
%
% INPUTS
%
%
% OUTPUTS
%
% u: list of physical commands of dim n [rad/s].

u = {0*pi/180,1.5*pi/180,-1.5*pi/180,3*pi/180,-3*pi/180};

end
