function s = vnncomp_status_to_str(code)
% vnncomp_status_to_str  Map an integer verification status code to its string.
%
%   Extracted (behavior-preserving) copy of the local `status_to_str` helper in
%   run_all_benchmarks.m so it can be unit-tested in isolation. Adopting it as the
%   single source of truth in run_all_benchmarks.m is a recommended follow-up
%   refactor (delete the local function and call this instead); until then this
%   file MUST be kept byte-for-byte equivalent to that local function.
%
%   Status codes:
%      0 = sat      1 = unsat       2 = unknown
%     -1 = error   -2 = missing    -3 = decompress_failed
%   anything else -> sprintf('code_%d', code)

switch code
    case 0,  s = 'sat';
    case 1,  s = 'unsat';
    case 2,  s = 'unknown';
    case -1, s = 'error';
    case -2, s = 'missing';
    case -3, s = 'decompress_failed';
    otherwise, s = sprintf('code_%d', code);
end
end
