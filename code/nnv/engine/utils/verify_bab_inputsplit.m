function [status, cex] = verify_bab_inputsplit(nnvnet, lb, ub, prop, opts)
% VERIFY_BAB_INPUTSPLIT  Input-domain branch-and-bound verification.
%
%   [status, cex] = verify_bab_inputsplit(nnvnet, lb, ub, prop, opts)
%
% Recursively bisects the input box and runs the EXISTING sound approx-star
% pipeline on each cell: cells proven safe are dropped; undecided cells are
% split along their longest normalized edge; each popped cell first gets a
% cheap concrete falsification probe (corners + center + random samples) so
% SAT exits early with a witness. This is T1 of REACH_PERFORMANCE_STRATEGIES.md
% -- the input-space BaB route the rest of the field (PyRAT, CORA,
% alpha-beta-CROWN's input-split mode) uses to crack low-input-dimension
% benchmarks like acasxu (5-D), where NNV's exact-star path enumeration pays
% ~143x more per path (CAV 2020, Table 1) and times out.
%
% Inputs:
%   nnvnet : NNV NN object (matlab2nnv output) -- reach + evaluate used
%   lb, ub : flat input bounds (vnnlib order), column vectors
%   prop   : cell array of specs, each a struct with .Hg = HalfSpace array
%            (the UNSAFE region; same contract as run_vnncomp_instance)
%   opts   : struct with fields
%       .timeout      total budget in seconds (driver returns 2 on expiry)
%       .makeSet      @(lb,ub) -> input set for nnvnet.reach (lets the caller
%                     apply its create_input_set reshape/ImageStar conventions)
%       .evalSample   @(x_flat) -> network output for a concrete flat input
%                     (caller applies the same reshape conventions)
%       .reachOptions reach options struct for the per-cell reach
%                     (default: approx-star)
%       .nSamples     random samples per cell in the falsification probe
%                     (default 5; corners+center are always tried, corners
%                     capped at 2^10)
%       .maxCells     safety cap on processed cells (default 20000)
%
% Outputs:
%   status : 1 = UNSAT proven (every cell sound-safe)
%            0 = SAT, cex = {x_flat; y} concrete validated witness
%            2 = unknown (timeout / cell cap / any reach failure on an
%                unsplittable cell)
%
% SOUNDNESS: the cells exactly tile the input box (bisection at midpoints,
% closed cells -- boundary overlap is sound for safety proofs). UNSAT is
% only returned when EVERY cell is proven safe by the sound approx-star
% reach + verify_specification. SAT is only returned for a CONCRETE input
% whose output violates the spec (the caller additionally replays witnesses
% through onnxruntime before emitting `sat`). Everything else is unknown.

    if ~isfield(opts, 'reachOptions') || isempty(opts.reachOptions)
        opts.reachOptions = struct('reachMethod', 'approx-star');
    end
    if ~isfield(opts, 'nSamples'), opts.nSamples = 5; end
    if ~isfield(opts, 'maxCells'), opts.maxCells = 20000; end

    lb = double(lb(:)); ub = double(ub(:));
    width0 = max(ub - lb, eps);   % normalization for the split heuristic

    % worklist of cells; process widest-first (max normalized edge) so the
    % hard regions split early and easy regions clear in large chunks
    Q = struct('lb', lb, 'ub', ub);
    t0 = tic; processed = 0;
    status = 1; cex = [];

    while ~isempty(Q)
        if toc(t0) > opts.timeout || processed >= opts.maxCells
            status = 2; cex = []; return;
        end
        % pop the cell with the largest normalized width
        [~, qi] = max(arrayfun(@(c) max((c.ub - c.lb) ./ width0), Q));
        c = Q(qi); Q(qi) = [];
        processed = processed + 1;

        % 1) cheap concrete falsification: corners + center + random samples
        w = quick_falsify(c.lb, c.ub, prop, opts);
        if ~isempty(w)
            status = 0; cex = w; return;
        end

        % 2) sound reach on the cell
        cellSafe = false;
        try
            IS = opts.makeSet(c.lb, c.ub);
            R = nnvnet.reach(IS, opts.reachOptions);
            res = verify_specification(R, prop);
            cellSafe = (res == 1);
        catch
            % reach failed on this cell: try splitting further (smaller cells
            % are often easier); if it cannot be split, give up soundly
        end

        if ~cellSafe
            [c1, c2, splittable] = bisect_longest_edge(c, width0);
            if ~splittable
                status = 2; cex = []; return;   % cannot refine -> unknown (sound)
            end
            Q(end+1) = c1; %#ok<AGROW>
            Q(end+1) = c2; %#ok<AGROW>
        end
    end
    % worklist empty: every cell proven safe -> UNSAT
end

% -------------------------------------------------------------------------

function w = quick_falsify(lb, ub, prop, opts)
    % concrete probe: center + corners (capped) + uniform random samples.
    % returns {x; y} on a genuine violation, [] otherwise.
    w = [];
    n = numel(lb);
    X = (lb + ub) / 2;                       % center
    if n <= 10                               % all corners only for small dims
        corners = dec2bin(0:2^n-1) - '0';
        X = [X, lb + corners' .* (ub - lb)];
    end
    if opts.nSamples > 0
        X = [X, lb + rand(n, opts.nSamples) .* (ub - lb)];
    end
    for i = 1:size(X, 2)
        x = X(:, i);
        try
            y = opts.evalSample(x);
        catch
            continue;
        end
        y = double(y(:));
        for p = 1:numel(prop)
            Hs = prop{p}.Hg;
            for d = 1:numel(Hs)
                if all(Hs(d).G * y <= Hs(d).g + 1e-12)
                    w = {x; y};              % concrete unsafe point
                    return;
                end
            end
        end
    end
end

function [c1, c2, ok] = bisect_longest_edge(c, width0)
    % split the cell at the midpoint of its longest NORMALIZED edge
    wn = (c.ub - c.lb) ./ width0;
    [wmax, k] = max(wn);
    mid = (c.lb(k) + c.ub(k)) / 2;
    ok = wmax > 0 && mid > c.lb(k) && mid < c.ub(k);   % representable progress
    c1 = c; c2 = c;
    if ok
        c1.ub(k) = mid;
        c2.lb(k) = mid;
    end
end
