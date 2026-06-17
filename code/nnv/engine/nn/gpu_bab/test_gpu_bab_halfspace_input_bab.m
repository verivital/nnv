function tests = test_gpu_bab_halfspace_input_bab
% TEST_GPU_BAB_HALFSPACE_INPUT_BAB  Soundness of the input-bisection halfspace BaB on a tiny
%   bilinear net (output y = x1*x2 over [-1,1]^2, so y in [-1,1]):
%     - a genuinely SAFE unsafe-region (y >= 2, empty over the box) must certify -> 'robust';
%     - a genuinely REACHABLE unsafe-region (y <= 0.5) must NOT be falsely certified -> not 'robust'.
%   This exercises the product McCormick relaxation + the bisection + the certify path, and guards
%   the soundness invariant (never a false 'robust'). Run: runtests('test_gpu_bab_halfspace_input_bab').
    tests = functiontests(localfunctions);
end

function test_safe_certifies_robust(tc)
    ops = i_net();
    lb = [-1; -1]; ub = [1; 1];
    Gd = {[-1]}; gd = {[-2]};      % unsafe disjunct {-y <= -2} == {y >= 2}: empty over the box -> SAFE
    [v, info] = gpu_bab_halfspace_input_bab(ops, lb, ub, Gd, gd, ...
        struct('maxNodes', 20000, 'minWidth', 1e-7, 'timeCap', 60));
    verifyEqual(tc, v, 'robust', sprintf('safe spec not certified (reason=%s)', info.reason));
end

function test_unsafe_not_falsely_robust(tc)
    ops = i_net();
    lb = [-1; -1]; ub = [1; 1];
    Gd = {[1]}; gd = {[0.5]};      % unsafe disjunct {y <= 0.5}: reachable (e.g. x=(-1,1)->y=-1) -> NOT safe
    [v, ~] = gpu_bab_halfspace_input_bab(ops, lb, ub, Gd, gd, ...
        struct('maxNodes', 20000, 'minWidth', 1e-3, 'timeCap', 60));
    verifyNotEqual(tc, v, 'robust', 'UNSOUND: falsely certified a reachable unsafe region as robust');
end

function ops = i_net()
% y = x1 * x2 : two affine selectors feeding an elementwise product.
    ops = { struct('type','affine','W',[1 0],'b',0,'src',0), ...   % op1 -> x1
            struct('type','affine','W',[0 1],'b',0,'src',0), ...   % op2 -> x2
            struct('type','product','inputs',[1 2],'sizes',[1 1],'src',1) };  % op3 -> x1*x2
end
