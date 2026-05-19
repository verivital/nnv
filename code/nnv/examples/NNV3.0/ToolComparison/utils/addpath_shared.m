function addpath_shared()
%ADDPATH_SHARED  Add the v2 utils dir to the MATLAB path.
%   After v2 was self-contained (tool_utils / rebuild_for_aivl /
%   toolbox_install / parse_argmax_vnnlib all copied locally), this is
%   just an addpath of the same directory the caller already lives in.
%   Kept as a function for back-compat with any code that still calls it.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
end
