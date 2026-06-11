function [status, elapsed, errmsg, out] = reach_with_timeout(net, IS, timeout_sec, reachOptions)
%REACH_WITH_TIMEOUT  Run net.reach with a hard wall-clock timeout via parfeval.
%   status:  'OK', 'TIMEOUT', 'FAIL'
%
% Uses an EXISTING parallel pool if one is up. Caller should ensure a pool
% exists (pp = gcp('nocreate'); or parpool('local', 1)) once at startup so
% individual reach calls don't pay the pool-startup cost.
%
% If no pool is available, falls back to a try/catch with a soft budget
% (cannot mid-call cancel) — the elapsed time is reported but the process
% is not killed.

if nargin < 3 || isempty(timeout_sec), timeout_sec = 60; end
if nargin < 4, reachOptions = struct('reachMethod', 'approx-star'); end

out = [];
errmsg = '';
t = tic;

pp = gcp('nocreate');
if ~isempty(pp)
    f = parfeval(pp, @(net_, IS_, opts_) net_.reach(IS_, opts_), 1, net, IS, reachOptions);
    ok_wait = wait(f, 'finished', timeout_sec);
    if ~ok_wait
        cancel(f);
        elapsed = timeout_sec;
        status = 'TIMEOUT';
        return;
    end
    if ~isempty(f.Error)
        elapsed = toc(t);
        status = 'FAIL';
        errmsg = f.Error.message;
        return;
    end
    out = fetchOutputs(f);
    elapsed = toc(t);
    status = 'OK';
    return;
end

% No pool: best-effort try/catch with soft budget
try
    out = net.reach(IS, reachOptions);
    elapsed = toc(t);
    status = 'OK';
catch ME
    elapsed = toc(t);
    status = 'FAIL';
    errmsg = ME.message;
end
end
