function tests = test_LstmLayer_evaluate_outputMode
%TEST_LSTMLAYER_EVALUATE_OUTPUTMODE  Regression tests for GitHub issues #325 and #326.
%
%   #325: LstmLayer.evaluate (and evaluateSequence, which delegates to it)
%   discarded the first output of lstm() and returned only the final hidden
%   state regardless of OutputMode. With OutputMode='sequence' they must
%   return the full hidden-state sequence (numHiddenUnits x seqLength),
%   matching MATLAB's lstmLayer semantics and the layer's own star
%   reachability path (reach_star_single_input, outputMode == "sequence",
%   fixed in PR #322). OutputMode='last' behavior is unchanged.
%
%   #326: LstmLayer.reach_zono was non-functional (referenced nonexistent
%   properties obj.Weights/obj.Bias/obj.OutputSize, wrong input-size
%   semantics) and a single affine map has no sound LSTM semantics, so the
%   'approx-zono' dispatch in reach()/reachSequence() must fail closed with
%   error id NNV:LstmLayer:zonoUnsupported (an error is safe; a wrong set
%   is catastrophic).
%
%   To run: results = runtests('test_LstmLayer_evaluate_outputMode')
    tests = functiontests(localfunctions);
end

% ---------- #325: evaluate/evaluateSequence OutputMode handling ----------

function test_sequence_mode_shape_and_last_column(tc)
    % OutputMode='sequence' must return numHiddenUnits x seqLength, and its
    % last column must equal the OutputMode='last' result on the same input
    [Lseq, Llast, inputSize, numHiddenUnits, seqLength] = make_lstm_pair(42);

    x = randn(inputSize, seqLength);
    yseq = Lseq.evaluateSequence(x);
    ylast = Llast.evaluateSequence(x);

    verifyEqual(tc, size(yseq), [numHiddenUnits, seqLength], ...
        "OutputMode='sequence' must return the full hidden-state sequence (features x time)");
    verifyEqual(tc, size(ylast), [numHiddenUnits, 1], ...
        "OutputMode='last' must return only the final hidden state");
    verifyEqual(tc, double(yseq(:, end)), double(ylast), 'AbsTol', 1e-12, ...
        "last column of the sequence must equal the OutputMode='last' result");

    % evaluate and evaluateSequence must agree (the latter delegates)
    verifyEqual(tc, double(Lseq.evaluate(x)), double(yseq), 'AbsTol', 1e-12, ...
        'evaluate and evaluateSequence must agree');
    verifyEqual(tc, double(Llast.evaluate(x)), double(ylast), 'AbsTol', 1e-12, ...
        'evaluate and evaluateSequence must agree');
end

function test_sequence_mode_prefix_consistency(tc)
    % column t of the sequence output equals the 'last' output on the
    % prefix x(:, 1:t) (the initial hidden/cell states are fixed), pinning
    % the values of every column, not just the final one
    [Lseq, Llast, inputSize, ~, seqLength] = make_lstm_pair(43);

    x = randn(inputSize, seqLength);
    yseq = Lseq.evaluateSequence(x);
    for t = 1:seqLength
        ht = Llast.evaluateSequence(x(:, 1:t));
        verifyEqual(tc, double(yseq(:, t)), double(ht), 'AbsTol', 1e-9, ...
            sprintf('column %d must equal the final hidden state of the prefix x(:,1:%d)', t, t));
    end
end

% ---------- #325: evaluate-vs-reach consistency (the core divergence) ----------

function test_reach_consistency_sequence_mode(tc)
    % for a degenerate (zero-width) ImageStar input, the center of the
    % reachSequence output must equal the concrete evaluation; this catches
    % evaluate-vs-reach divergence (pre-fix, evaluateSequence returned
    % numHiddenUnits x 1 while the reach path returned one column per step)
    [Lseq, ~, inputSize, numHiddenUnits, seqLength] = make_lstm_pair(44);

    x = randn(inputSize, seqLength);
    R = Lseq.reachSequence(ImageStar(x, x), 'approx-star');
    verifyClass(tc, R, 'ImageStar');
    verifyEqual(tc, R.height, numHiddenUnits, 'output must have numHiddenUnits rows');
    verifyEqual(tc, R.width, seqLength, ...
        'sequence-mode output must have one column per time step');

    [lb, ub] = R.getRanges;
    c = reshape((lb + ub)/2, numHiddenUnits, seqLength);
    y = Lseq.evaluateSequence(x);
    verifyEqual(tc, size(y), size(c), ...
        'evaluateSequence output shape must match the reach output shape');
    verifyEqual(tc, double(c), double(y), 'AbsTol', 1e-6, ...
        'reach center must equal the concrete evaluation for a point input');
    verifyLessThan(tc, max(abs(ub(:) - lb(:))), 1e-6, ...
        'reach output for a point input must be (numerically) a point');
end

function test_reach_consistency_last_mode(tc)
    % same degenerate-input consistency for OutputMode='last' (guards the
    % byte-identical 'last' behavior against regressions)
    [~, Llast, inputSize, numHiddenUnits, seqLength] = make_lstm_pair(45);

    x = randn(inputSize, seqLength);
    R = Llast.reachSequence(ImageStar(x, x), 'approx-star');
    verifyClass(tc, R, 'ImageStar');
    verifyEqual(tc, R.height, numHiddenUnits, 'output must have numHiddenUnits rows');

    [lb, ub] = R.getRanges;
    c = reshape((lb + ub)/2, numHiddenUnits, 1);
    y = Llast.evaluateSequence(x);
    verifyEqual(tc, size(y), size(c), ...
        'evaluateSequence output shape must match the reach output shape');
    verifyEqual(tc, double(c), double(y), 'AbsTol', 1e-6, ...
        'reach center must equal the concrete evaluation for a point input');
end

% ---------- #326: approx-zono must fail closed ----------

function test_reach_zono_fails_closed(tc)
    % the zonotope path is non-functional and has no sound semantics; both
    % dispatch entry points must raise a clear, identifiable error instead
    % of returning a wrong set (VNN-COMP: an error is safe, a wrong set is
    % a -150 catastrophe)
    [Lseq, Llast, inputSize, ~, seqLength] = make_lstm_pair(46);

    x = randn(inputSize, seqLength);
    IZ = ImageZono(x - 0.1, x + 0.1);

    verifyError(tc, @() Llast.reach(IZ, 'approx-zono', []), ...
        'NNV:LstmLayer:zonoUnsupported');
    verifyError(tc, @() Llast.reachSequence(IZ, 'approx-zono', []), ...
        'NNV:LstmLayer:zonoUnsupported');
    verifyError(tc, @() Lseq.reach(IZ, 'approx-zono', []), ...
        'NNV:LstmLayer:zonoUnsupported');
    verifyError(tc, @() Lseq.reachSequence(IZ, 'approx-zono', []), ...
        'NNV:LstmLayer:zonoUnsupported');

    % direct calls to the stubbed methods fail closed as well
    verifyError(tc, @() Llast.reach_zono(IZ), ...
        'NNV:LstmLayer:zonoUnsupported');
    verifyError(tc, @() Llast.reach_zono_multipleInputs(IZ, []), ...
        'NNV:LstmLayer:zonoUnsupported');
end

% ---------- helpers ----------

function [Lseq, Llast, inputSize, numHiddenUnits, seqLength] = make_lstm_pair(seed)
    % two small deterministic LstmLayers sharing identical fixed-rng
    % weights, differing only in OutputMode (constructor pattern mirrors
    % test_soundness_LstmLayer_reachSequence.m)
    rng(seed);

    inputSize = 3;
    numHiddenUnits = 4;
    seqLength = 5;

    inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
    recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
    bias = 0.1*randn(4*numHiddenUnits, 1);
    cellState = zeros(numHiddenUnits, 1);
    hiddenState = zeros(numHiddenUnits, 1);

    Lseq = LstmLayer('Name', 'lstm_eval_seq', ...
        'InputSize', inputSize, ...
        'NumHiddenUnits', numHiddenUnits, ...
        'InputWeights', inputWeights, ...
        'RecurrentWeights', recurrentWeights, ...
        'Bias', bias, ...
        'CellState', cellState, ...
        'HiddenState', hiddenState, ...
        'OutputMode', 'sequence');

    Llast = LstmLayer('Name', 'lstm_eval_last', ...
        'InputSize', inputSize, ...
        'NumHiddenUnits', numHiddenUnits, ...
        'InputWeights', inputWeights, ...
        'RecurrentWeights', recurrentWeights, ...
        'Bias', bias, ...
        'CellState', cellState, ...
        'HiddenState', hiddenState, ...
        'OutputMode', 'last');
end
