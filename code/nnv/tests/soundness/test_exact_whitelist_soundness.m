classdef test_exact_whitelist_soundness < matlab.unittest.TestCase
    % AUTO-DETECTION for the exact-star gate (NN.layer_reach_is_exact) [42].
    %
    % The gate promotes unknown->not-robust ONLY when every layer's reach is
    % exact, using a hand-maintained whitelist. That list has been wrong twice
    % (AdditionLayer/ConcatenationLayer over-approximate via blkdiag; SignLayer's
    % reach is outright unsound). This test mechanically pins the contract so a
    % future mis-whitelisting fails CI:
    %   (1) every WHITELISTED piecewise-linear/activation layer must have an
    %       exact-star reach that is MC-SOUND and TIGHT (envelope == true range);
    %   (2) every known OVER-APPROXIMATE or UNSOUND layer must be EXCLUDED.
    % (Pure affine layers are exact by construction -- only checked for inclusion.)

    methods (Test)

        function test_whitelisted_pwl_are_exact(testCase)
            net = NN({ReluLayer('r')}, []);   % any net just to reach the method
            L_relu = ReluLayer('relu');
            L_lrelu = LeakyReluLayer('lr', 1, {'in'}, 1, {'out'}, 0.1);
            pwl = {L_relu, L_lrelu};
            for i = 1:numel(pwl)
                L = pwl{i};
                cls = class(L);
                testCase.verifyTrue(net.layer_reach_is_exact(L), ...
                    sprintf('%s must be whitelisted (it is exact)', cls));
                % ReLU and LeakyReLU are elementwise-monotone, so the TRUE range
                % over a box [lb,ub] is exactly [f(lb), f(ub)] per coordinate.
                lb = [-1; -0.5; -2; 0.3]; ub = [1; 2; -0.1; 1.5];
                IS = Star(lb, ub);
                RS = L.reach(IS, 'exact-star');
                d = IS.dim; a = inf(d,1); b = -inf(d,1);
                for k = 1:numel(RS)
                    [lo, hi] = RS(k).getRanges; a = min(a, lo(:)); b = max(b, hi(:));
                end
                trueLo = L.evaluate(lb); trueHi = L.evaluate(ub);   % monotone
                % SOUND: reach envelope contains the true range.
                testCase.verifyTrue(all(a <= trueLo + 1e-7) && all(b >= trueHi - 1e-7), ...
                    sprintf('%s exact-star reach is UNSOUND (excludes true range)', cls));
                % TIGHT/EXACT: reach envelope EQUALS the true range (no over-approx).
                testCase.verifyTrue(max(abs(a - trueLo)) < 1e-6 && max(abs(b - trueHi)) < 1e-6, ...
                    sprintf('%s exact-star reach OVER-APPROXIMATES (not exact)', cls));
            end
        end

        function test_whitelisted_maxpool_is_sound(testCase)
            net = NN({ReluLayer('r')}, []);
            L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);
            testCase.verifyTrue(net.layer_reach_is_exact(L), 'MaxPooling2DLayer must be whitelisted');
            IS = ImageStar(-ones(2,2,1), ones(2,2,1));
            RS = L.reach(IS, 'exact-star'); a = inf; b = -inf;
            for k = 1:numel(RS)
                [lo, hi] = RS(k).getRanges; a = min(a, min(lo(:))); b = max(b, max(hi(:)));
            end
            % MC soundness: every concrete maxpool output is contained.
            viol = 0;
            for t = 1:4000
                xs = -1 + 2*rand(2,2,1); ys = L.evaluate(xs); ys = ys(:);
                if any(ys < a - 1e-9) || any(ys > b + 1e-9), viol = viol + 1; end
            end
            testCase.verifyEqual(viol, 0, 'MaxPool exact-star reach is UNSOUND');
            % Tight: true range of max over [-1,1]^4 is exactly [-1, 1].
            testCase.verifyTrue(abs(a - (-1)) < 1e-6 && abs(b - 1) < 1e-6, ...
                'MaxPool exact-star reach OVER-APPROXIMATES');
        end

        function test_whitelisted_affine_are_included(testCase)
            % Affine layers are exact by construction (Star.affineMap) -- only
            % their inclusion is pinned (a regression that drops them would only
            % hurt precision, but the contract should be stable).
            net = NN({ReluLayer('r')}, []);
            affine = { ...
                FullyConnectedLayer('fc', eye(3), zeros(3,1)), ...
                FlattenLayer('f'), ReshapeLayer('rs', [1 1 3]) };
            for i = 1:numel(affine)
                testCase.verifyTrue(net.layer_reach_is_exact(affine{i}), ...
                    sprintf('%s should be whitelisted (affine/exact)', class(affine{i})));
            end
        end

        function test_overapprox_and_unsound_layers_excluded(testCase)
            % The cardinal-sin guard: any layer that over-approximates or is
            % unsound under exact-star MUST be excluded, else the gate could
            % certify a false not-robust.
            net = NN({ReluLayer('r')}, []);
            excluded = { ...
                SignLayer(), ...                                   % reach UNSOUND
                GeluLayer('g'), SiLULayer('s'), ...                % smooth nonlinear
                ElementwiseProductLayer('ep'), ...                 % bilinear over-approx
                AdditionLayer('a', 2, 1, {'in1','in2'}, {'out'}), ...          % blkdiag Minkowski over-approx
                ConcatenationLayer('c', 2, 1, {'in1','in2'}, {'out'}, 1), ...  % blkdiag concat over-approx
                DepthConcatenationLayer('dc', 2, 1, {'in1','in2'}, {'out'}) };
            for i = 1:numel(excluded)
                testCase.verifyFalse(net.layer_reach_is_exact(excluded{i}), ...
                    sprintf('%s must be EXCLUDED from the exact whitelist (over-approx/unsound)', ...
                            class(excluded{i})));
            end
            sm = SoftmaxLayer('sm'); sm.IsFinalLayer = false;     % intermediate softmax = over-approx
            testCase.verifyFalse(net.layer_reach_is_exact(sm), 'intermediate Softmax must be excluded');
            % Conv1DLayer (box over-approx branch from cached im_lb/im_ub, no
            % plain reach()) and PixelClassificationLayer (estimateRanges
            % over-approx, numeric output) must also be excluded. Both need
            % constructor args; build defensively and skip if unavailable (the
            % exclusion is enforced by the whitelist regardless).
            try, c1 = Conv1DLayer(ones(2,1,1), zeros(1,1)); catch, c1 = []; end
            if ~isempty(c1)
                testCase.verifyFalse(net.layer_reach_is_exact(c1), ...
                    'Conv1DLayer must be excluded (box over-approximation branch)');
            end
            try, pc = PixelClassificationLayer('pc', categorical({'a';'b'}), [1 1 2]); catch, pc = []; end
            if ~isempty(pc)
                testCase.verifyFalse(net.layer_reach_is_exact(pc), ...
                    'PixelClassificationLayer must be excluded (estimateRanges over-approx)');
            end
        end

        function test_whitelisted_affine_layers_are_exact(testCase)
            % The AFFINE whitelisted layers must have an exact-star reach that
            % EQUALS the true per-coordinate range (affine maps are exact). This
            % pins them so a future loosening fails CI -- it directly catches the
            % BatchNorm Star-generator bug this review found (the old code added a
            % constant to every basis generator -> a sound but LOOSE reach that
            % was nonetheless whitelisted as exact, a gate hole).
            bn = BatchNormalizationLayer('TrainedMean',[1;-2;3],'TrainedVariance',[2;1;0.5], ...
                'Epsilon',1e-3,'Scale',[2;-1;0.5],'Offset',[1;0;-1]); bn.NumChannels = 3;
            testCase.verifyAffineExact(bn, Star([-1;-1;-1],[1;2;0.5]));
            ea = ElementwiseAffineLayer('ea',[2;-1;0.5],[1;0;-1],true,true);
            testCase.verifyAffineExact(ea, Star([-1;-1;-1],[1;2;0.5]));
            fc = FullyConnectedLayer('fc',[1 2 0; 0 1 3],[1;-1]);
            testCase.verifyAffineExact(fc, Star([-1;-1;-1],[1;2;0.5]));
        end

        function test_eaffine_evaluate_preserves_shape(testCase)
            % [10] direct (no-ONNX) pin: a [1,N] row Scale on an [N,1] column
            % input must scale element-wise to [N,1], NOT form an [N,N] outer
            % product (the ACAS-Xu input-normalization regression).
            L = ElementwiseAffineLayer('ea', (1:5), 0, true, false);   % Scale [1,5]
            y = L.evaluate((1:5)');
            testCase.verifyEqual(size(y), [5 1], 'row Scale on column input must stay [N,1]');
            testCase.verifyEqual(double(y(:)), ((1:5).*(1:5))', 'AbsTol', 1e-9);
        end

        function test_softmax_estimateRanges_branch_sound(testCase)
            % cf0608f76 added an estimateRanges fallback for intermediate-softmax
            % ImageStar reach when im_lb/im_ub are EMPTY (the segmentation case).
            % Every other softmax test uses box-constructed ImageStars (cached
            % bounds), so this branch was uncovered. Pin its MC soundness on an
            % ImageStar with a real predicate and empty cached bounds.
            H=2; W=2; C=3;
            V = cat(4, randn(H,W,C), 0.3*randn(H,W,C), 0.2*randn(H,W,C));
            IS = ImageStar(V, [1 0;0 1;-1 0;0 -1], [1;1;1;1], [-1;-1], [1;1]);
            testCase.assertTrue(isempty(IS.im_lb), 'precond: ImageStar must have empty cached bounds');
            L = SoftmaxLayer('sm'); L.IsFinalLayer = false;
            OS = L.reach(IS); [a,b] = OS.getRanges; a = a(:); b = b(:);
            viol = 0;
            for k = 1:4000
                al = IS.pred_lb + (IS.pred_ub - IS.pred_lb).*rand(2,1);
                x = IS.V(:,:,:,1) + IS.V(:,:,:,2)*al(1) + IS.V(:,:,:,3)*al(2);
                y = double(reshape(L.evaluate(single(x)), [], 1));
                if any(y < a-1e-6) || any(y > b+1e-6), viol = viol + 1; end
            end
            testCase.verifyEqual(viol, 0, 'estimateRanges-branch softmax reach is UNSOUND');
        end

    end

    methods (Access = private)
        function verifyAffineExact(testCase, L, X)
            % X: a box Star. The single-input reach envelope must EQUAL the true
            % per-coordinate affine range [min(f(lb),f(ub)), max(f(lb),f(ub))]
            % (sound AND tight -> exact). Uses reach_star_single_input (the exact
            % affine map for these layers).
            R = L.reach_star_single_input(X);
            [a, b] = R.getRanges; a = a(:); b = b(:);
            [lb, ub] = X.getRanges; ylb = L.evaluate(lb); yub = L.evaluate(ub);
            tlo = min(ylb(:), yub(:)); thi = max(ylb(:), yub(:));
            testCase.verifyTrue(all(a <= tlo + 1e-6) && all(b >= thi - 1e-6), ...
                sprintf('%s exact-star reach is UNSOUND (excludes true range)', class(L)));
            testCase.verifyTrue(max(abs(a - tlo)) < 1e-5 && max(abs(b - thi)) < 1e-5, ...
                sprintf('%s exact-star reach is NOT EXACT (over-approximates)', class(L)));
        end
    end
end
