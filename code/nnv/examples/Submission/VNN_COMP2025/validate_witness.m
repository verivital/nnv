function ok = validate_witness(net, x, lb, ub, Hs, inputSize, inputFormat, in_tol, out_tol)
%VALIDATE_WITNESS  Confirm a SAT counterexample is REAL before emitting `sat`.
%   VNN-COMP scores an incorrect/non-replayable verdict at -150 (16x the +10 for
%   a correct one), and NNV leaked 19 such "incorrect/missing-CE" points in 2025.
%   This re-evaluates a candidate counterexample on the network and confirms it
%   genuinely (1) respects the input box and (2) lands in the unsafe output
%   region, so the caller never writes `sat` for a numerical near-miss. If it
%   fails, the caller must emit `unknown` instead. See VNNCOMP2026_STRATEGY.md
%   (Pillar 2).
%
%   Inputs:
%     net          dlnetwork OR an NNV NN (the same object falsification used).
%     x            candidate input, FLAT column vector (same order as lb/ub).
%     lb, ub       input box bounds (flat vectors).
%     Hs           array of HalfSpace (unsafe region); a real CE satisfies some h.
%     inputSize    network input size for reshaping x.
%     inputFormat  'default' or a dlarray dims label (for dlnetwork forward).
%     in_tol       input-constraint slack (default 1e-6; rules require ~zero tol).
%     out_tol      output-constraint slack (default 1e-4 absolute).
%
%   Output: ok (logical) -- true only if the witness is a validated counterexample.

    if nargin < 8 || isempty(in_tol),  in_tol  = 1e-6; end
    if nargin < 9 || isempty(out_tol), out_tol = 1e-4; end
    ok = false;

    x  = double(x(:));
    lb = double(lb(:)); ub = double(ub(:));

    % (1) Input must lie within the specified box (near-zero tolerance).
    if numel(x) ~= numel(lb) || any(x < lb - in_tol) || any(x > ub + in_tol)
        return;
    end

    % (2) Re-evaluate the network on the witness.
    try
        y = eval_net(net, x, inputSize, inputFormat);
    catch
        return;   % cannot validate -> treat as invalid (caller -> unknown)
    end
    y = double(y(:));
    if isempty(y) || any(~isfinite(y))
        return;
    end

    % (3) Confirm the output lands in some unsafe HalfSpace (G*y <= g).
    for h = 1:numel(Hs)
        m = Hs(h).G * y - Hs(h).g(:);
        if all(m <= out_tol)
            ok = true; return;
        end
    end
end

% ---- helpers ----

function y = eval_net(net, x, inputSize, inputFormat)
    if isa(net, 'NN')                 % NNV native net (Python-importer manifest)
        y = net.evaluate(x);
        return;
    end
    % dlnetwork: reshape + label + predict.
    if isscalar(inputSize)
        xr = reshape(x, [inputSize, 1]);
    else
        xr = reshape(x, inputSize);
    end
    fmt = resolve_format(inputFormat, inputSize);
    y = extractdata(predict(net, dlarray(single(xr), fmt)));
end

function fmt = resolve_format(inputFormat, inputSize)
    if nargin >= 1 && ~isempty(inputFormat) && ~strcmp(inputFormat, "default")
        fmt = char(inputFormat);
        if ~endsWith(fmt, 'B'), fmt = [fmt 'B']; end
        return;
    end
    if isscalar(inputSize) || (numel(inputSize) <= 3 && nnz(inputSize > 1) <= 1)
        fmt = 'CB';
    else
        fmt = 'SSCB';
    end
end
