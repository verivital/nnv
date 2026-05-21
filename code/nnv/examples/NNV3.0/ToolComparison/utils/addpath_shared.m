function addpath_shared()
%ADDPATH_SHARED  Add the utils dir to the MATLAB path.
%   Self-contained (tool_utils / rebuild_for_aivl / parse_argmax_vnnlib
%   all live locally), so this is just an addpath of the same directory
%   the caller already lives in. Kept as a function for back-compat with
%   any code that still calls it.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
end
