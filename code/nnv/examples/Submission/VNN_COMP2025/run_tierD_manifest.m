function run_tierD_manifest(bench_root)
% RUN_TIERD_MANIFEST  Smoke-test the Tier-D benchmarks via the Python-importer
% MANIFEST path (load_nnv_from_mat) instead of the matlab2nnv category path.
%
% cgan / soundnessbench / test fail to import on the category (matlab2nnv) path
% (custom ReshapeLayer1000 / no dispatcher), but their ONNX is faithfully imported
% by tools/onnx2nnv_python/onnx2nnv.py into <model>.nnv.mat. This confirms that,
% loaded from the manifest, they forward-pass and produce a (sound) reach verdict.
%
%   run_tierD_manifest                       % default benchmark root
%   run_tierD_manifest('C:\path\to\benchmarks')

    if nargin < 1 || isempty(bench_root)
        script_dir = fileparts(mfilename('fullpath'));
        bench_root = fullfile(script_dir, '..','..','..','..','..','..', ...
            'vnncomp2025_benchmarks','benchmarks');
    end
    bench_root = char(bench_root);

    cases = {
        'cgan_2023',      'cGAN_imgSz32_nCh_1.onnx', 'cGAN_imgSz32_nCh_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib';
        'soundnessbench', 'model.onnx',              'model_0.vnnlib';
        'test',           'test_nano.onnx',          'test_nano.vnnlib';
    };

    fprintf('\n=== Tier-D manifest-path smoke test ===\n');
    for i = 1:size(cases,1)
        sub  = cases{i,1};
        onnx = fullfile(bench_root, sub, 'onnx',   cases{i,2});
        vnnl = fullfile(bench_root, sub, 'vnnlib', cases{i,3});
        manifest = regexprep(onnx, '\.onnx$', '.nnv.mat');
        fprintf('\n[%d] %s\n  manifest: %s\n', i, sub, manifest);
        if ~isfile(manifest)
            fprintf('  MISSING manifest (generate: python onnx2nnv.py %s)\n', onnx);
            continue;
        end
        try
            net = load_nnv_from_mat(manifest);
            property = load_vnnlib(vnnl);
            lb = property.lb; ub = property.ub; prop = property.prop;
            if iscell(lb), lb0 = lb{1}; ub0 = ub{1}; else, lb0 = lb; ub0 = ub; end
            p0 = prop{1};

            % forward-pass sanity at the box center
            yc = net.evaluate(double((lb0(:)+ub0(:))/2));
            fprintf('  forward OK (out dim %d)\n', numel(yc));

            IS = Star(double(lb0(:)), double(ub0(:)));
            st = 2;
            % Smoke test: a single COMPLETED reach confirms the manifest path works;
            % stop at the first method that returns (don't grind slower fallbacks).
            % approx-zono/abs-dom only matter if approx-star errors on this set type.
            for m = {'approx-star','abs-dom'}
                try
                    R  = net.reach(IS, struct('reachMethod', m{1}));
                    st = verify_specification(R, {p0});
                    fprintf('  reach %-12s -> %s\n', m{1}, statusStr(st));
                    break;   % completed -> done
                catch ME
                    fprintf('  reach %-12s errored: %s\n', m{1}, ME.message);
                end
            end
            fprintf('  RESULT (%s): %s\n', sub, statusStr(st));
        catch ME
            fprintf('  MANIFEST PATH ERROR: %s\n', ME.message);
        end
    end
    fprintf('\n=== Tier-D smoke test done ===\n');
end

function s = statusStr(st)
    switch st
        case 0, s = 'sat';
        case 1, s = 'unsat (SAFE)';
        case 2, s = 'unknown';
        otherwise, s = sprintf('code_%d', st);
    end
end
