classdef test_concat_mixed_empty_C < matlab.unittest.TestCase
    % Regression for the ConcatenationLayer mixed empty-C / non-empty-C bug
    % (GitHub Copilot review of PR #290).
    %
    % When one input to a Concat is a BOX (empty constraint matrix C) and another
    % is CONSTRAINED (non-empty C, e.g. from an approx-star ReLU reach), the old
    % code combined constraints with `blkdiag` over only the non-empty C's. If an
    % EARLIER input had empty C, its predicate columns were skipped, so new_C was
    % too narrow / mis-aligned vs the (correctly padded) new_V -> the Star/ImageStar
    % constructor rejected it ("Inconsistency between basic matrix and constraint
    % matrix"). Net effect: a perfectly valid box+constrained concatenation ERRORED
    % (fail-loud -> sound but a lost verdict) instead of producing the set.
    %
    % The fix pads each input's C to totalVars columns at its own predicate block.
    % This pins that the concat now (1) SUCCEEDS, (2) yields a totalVars-wide C, and
    % (3) is MC-SOUND (the reach set contains every concrete concatenation), for both
    % the Star (reach_concat_star) and ImageStar (reach_single_input) paths, and for
    % both input orderings (empty-C first AND second).

    methods (Test)

        function test_star_mixed_empty_C_sound(testCase)
            rng(0);
            for order = 1:2
                in1 = Star([0 1 0; 0 0 1], zeros(0,2), zeros(0,1), [-1;-1], [1;1]); % EMPTY C
                in2 = ReluLayer().reach(Star([-2;-2],[2;2]), 'approx-star');        % non-empty C
                if order == 1, ins = {in1, in2}; else, ins = {in2, in1}; end
                L = ConcatenationLayer('c', 2, 1, {'a','b'}, {'o'}, 1);

                R = L.reach_single_input(ins);   % must NOT error
                testCase.verifyEqual(size(R.C,2), R.nVar, ...
                    'combined C must be totalVars-wide');

                viol = testCase.mcViolationsStar(R, ins{1}, ins{2});
                testCase.verifyEqual(viol, 0, ...
                    sprintf('Star concat (order %d) reach is UNSOUND', order));
            end
        end

        function test_imagestar_mixed_empty_C_sound(testCase)
            rng(1);
            V = zeros(1,1,2,3); V(1,1,1,2) = 1; V(1,1,2,3) = 1;     % 2-channel box basis
            is1 = ImageStar(V, zeros(0,2), zeros(0,1), [-1;-1], [1;1]);          % EMPTY C
            is2 = ReluLayer().reach(ImageStar(-2*ones(1,1,2), 2*ones(1,1,2)), 'approx-star'); % non-empty C
            L = ConcatenationLayer('c', 2, 1, {'a','b'}, {'o'}, 3);  % concat on channel axis

            R = L.reach_single_input({is1, is2});   % must NOT error
            testCase.verifyEqual(size(R.C,2), R.numPred, ...
                'combined C must be totalVars-wide');

            viol = testCase.mcViolationsImageStar(R, is1, is2);
            testCase.verifyEqual(viol, 0, 'ImageStar concat reach is UNSOUND');
        end

    end

    methods (Access = private)

        function viol = mcViolationsStar(~, R, A, B)
            [lo,hi] = R.getRanges; lo = lo(:); hi = hi(:);
            viol = 0; got = 0; tries = 0;
            while got < 4000 && tries < 300000
                tries = tries + 1;
                aA = A.predicate_lb + (A.predicate_ub-A.predicate_lb).*rand(A.nVar,1);
                if ~isempty(A.C) && any(A.C*aA > A.d + 1e-9), continue; end
                aB = B.predicate_lb + (B.predicate_ub-B.predicate_lb).*rand(B.nVar,1);
                if ~isempty(B.C) && any(B.C*aB > B.d + 1e-9), continue; end
                y = [A.V(:,1)+A.V(:,2:end)*aA ; B.V(:,1)+B.V(:,2:end)*aB];
                got = got + 1;
                if any(y < lo-1e-6) || any(y > hi+1e-6), viol = viol + 1; end
            end
        end

        function viol = mcViolationsImageStar(~, R, A, B)
            [lo,hi] = R.getRanges; lo = lo(:); hi = hi(:);
            viol = 0; got = 0; tries = 0;
            while got < 4000 && tries < 300000
                tries = tries + 1;
                aA = A.pred_lb + (A.pred_ub-A.pred_lb).*rand(A.numPred,1);
                if ~isempty(A.C) && any(A.C*aA > A.d + 1e-9), continue; end
                aB = B.pred_lb + (B.pred_ub-B.pred_lb).*rand(B.numPred,1);
                if ~isempty(B.C) && any(B.C*aB > B.d + 1e-9), continue; end
                yA = A.V(:,:,:,1); for j = 1:A.numPred, yA = yA + A.V(:,:,:,j+1)*aA(j); end
                yB = B.V(:,:,:,1); for j = 1:B.numPred, yB = yB + B.V(:,:,:,j+1)*aB(j); end
                y = [yA(:); yB(:)];
                got = got + 1;
                if any(y < lo-1e-6) || any(y > hi+1e-6), viol = viol + 1; end
            end
        end

    end
end
