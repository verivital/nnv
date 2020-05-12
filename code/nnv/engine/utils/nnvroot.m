function [ nnvpath ] = nnvroot()
%NNVROOT
%   Returns the NNV root path (of the repository, so the repository /.)

    s = which('nnvroot');

    c = regexp(s, '[\\/]');
    if (size(c) < 2)
        error('path not found');
    end
    
    % this file should be located at code\nnv\engine\utils, so it is 4
    % levels above the git root
    fromBase = 4;

    nnvpath = s(1:c(end - fromBase) - 1);
end