%% Generic per-layer SOUNDNESS HARNESS (property-based MC-containment)
%
% Future-detection mechanism for reach() soundness regressions. For each
% registered single-input layer it draws ADVERSARIAL random input boxes
% (wide, nonuniform centers/radii, sign-straddling) plus corners+center, runs
% reach(), and asserts EVERY sampled evaluate() output lies in the reach set.
%
% Why this exists: the prior per-layer tests used tiny uniform boxes and missed
% real unsoundness (LayerNorm excluded reachable outputs by ~0.875; SiLU's error
% bound was sampling-only). Adversarial boxes expose those. A reach that is
% UNSOUND *or* that errors on a supported input fails this suite.
%
% Dev rule: adding/fixing any single-input reach() => register it here. This
% suite GATES CI (must never be allow-listed).

rng(7);
TOL = 1e-6;
NB = 60;    % adversarial boxes per layer
NS = 150;   % samples per box (corners + center + random)

% Each row: {name, layerObj, dim, centerScale, radiusMax}
% reach = layerObj.reach(Star,'approx-star'); oracle = layerObj.evaluate(x).
reg = {};
reg(end+1,:) = {'SiLULayer',    SiLULayer('silu'),                 5, 4, 6};
reg(end+1,:) = {'GeluLayer',    GeluLayer('gelu'),                 5, 4, 6};
reg(end+1,:) = {'SigmoidLayer', SigmoidLayer('sig'),               5, 4, 6};
reg(end+1,:) = {'TanhLayer',    TanhLayer('tanh'),                 5, 4, 6};
sm = SoftmaxLayer('sm'); sm.IsFinalLayer = false;
reg(end+1,:) = {'SoftmaxLayer(mid)', sm,                           5, 3, 5};
ln = LayerNormalizationLayer('Name','ln','NumChannels',5,'Epsilon',1e-5, ...
        'Scale',ones(5,1),'Offset',zeros(5,1));
reg(end+1,:) = {'LayerNormalizationLayer', ln,                     5, 3, 5};

fails = {};
for i = 1:size(reg,1)
    nm = reg{i,1}; L = reg{i,2}; d = reg{i,3}; cs = reg{i,4}; rmax = reg{i,5};
    [snd, nv, worst, wbox] = check_layer(L, d, cs, rmax, NB, NS, TOL);
    if snd
        fprintf('%-26s SOUND  (%d boxes x %d samples)\n', nm, NB, NS);
    else
        fprintf('%-26s UNSOUND: %d viol, worst overflow %.4g on box %s\n', nm, nv, worst, mat2str(wbox,3));
        fails{end+1} = sprintf('%s (%d viol, worst %.4g)', nm, nv, worst); %#ok<AGROW>
    end
end
assert(isempty(fails), 'Unsound layer reach detected:\n  %s', strjoin(fails, '\n  '));
fprintf('\n=== All registered layer reach()s are MC-sound on adversarial boxes ===\n');


function [snd, nv, worst, wbox] = check_layer(L, dim, cs, rmax, NB, NS, TOL)
    nv = 0; worst = 0; wbox = [];
    for t = 1:NB
        c = randn(dim,1)*cs;
        r = rand(dim,1)*rmax + 0.05;          % nonuniform radii
        lb = c - r; ub = c + r;
        try
            OS = L.reach(Star(lb, ub), 'approx-star');
        catch ME
            % A supported single-input reach that ERRORS is a harness failure
            % too (fail-loud is fine for unsupported ops, but these are supported).
            nv = nv + 1; worst = inf; wbox = [lb ub];
            fprintf('   (%s reach errored: %s)\n', class(L), ME.message);
            return;
        end
        [a, b] = OS.getRanges; a = a(:); b = b(:);
        for k = 1:NS
            if     k == 1, x = lb;
            elseif k == 2, x = ub;
            elseif k == 3, x = c;
            else,          x = lb + (ub - lb).*rand(dim,1);
            end
            y = L.evaluate(x); y = y(:);
            o = max(max(a - y), max(y - b));
            if o > TOL
                nv = nv + 1;
                if o > worst, worst = o; wbox = [lb ub]; end
            end
        end
    end
    snd = (nv == 0);
end
