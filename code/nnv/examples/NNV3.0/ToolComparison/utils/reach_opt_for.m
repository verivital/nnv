function reachOpt = reach_opt_for(alg)
%REACH_OPT_FOR  Map an algorithm string to an NNV reachOpt struct.
%   Caps exact-star parallelism at min(feature('numcores'), 8) — the default
%   Processes-profile pool size, and matches the inherited parfeval pool
%   used by the per-instance run_with_timeout wrapper (preventing
%   NN.start_pool from trying to delete+recreate the pool from a worker).
%
%   Algorithm strings:
%     approx-star
%     exact-star
%     relax-star-range-50, relax-star-range-75   (FC convention)
%     relax-star-area-50,  relax-star-area-75    (image convention)

    reachOpt = struct;
    reachOpt.reachOption = "single";
    reachOpt.numCores    = 1;
    reachOpt.relaxFactor = 0;
    switch alg
        case 'approx-star'
            reachOpt.reachMethod = 'approx-star';
        case 'relax-star-range-50'
            reachOpt.reachMethod = 'relax-star-range';
            reachOpt.relaxFactor = 0.5;
        case 'relax-star-range-75'
            reachOpt.reachMethod = 'relax-star-range';
            reachOpt.relaxFactor = 0.75;
        case 'relax-star-area-50'
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 0.5;
        case 'relax-star-area-75'
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 0.75;
        case 'exact-star'
            reachOpt.reachMethod = 'exact-star';
            reachOpt.reachOption = "parallel";
            % 8 matches the default Processes-profile pool size. The driver
            % wraps each call in parfeval on the client's pool; if exact-star
            % asks for a different size, NN.start_pool tries to delete+recreate
            % the pool from inside the worker and fails ("parpool cannot be
            % started from a worker"). 8-wide keeps the inherited pool reusable.
            reachOpt.numCores    = min(feature('numcores'), 8);
        otherwise
            error('ToolComparison:reach_opt_for:unknown_alg', ...
                  'Unknown NNV algorithm: %s', alg);
    end
end
