% test_load_vnnlib_forms
% Regression tests pinning the SHAPE of the property structure load_vnnlib
% produces for the common VNN-LIB output-specification forms used by the
% VNN-COMP benchmarks. These complement the single regression in
% test_load_vnnlib.m ("a constraint after an (or ...) is CONJOINED ...").
%
% Output-structure contract (see load_vnnlib.m header + HalfSpace.m):
%   property.prop                -> 1xN cell; each prop{n} is a struct with field .Hg
%   property.prop{n}.Hg          -> array of HalfSpace. An ARRAY (numel>1) is an
%                                   OR / disjunction. Within ONE HalfSpace, the rows
%                                   of .G (k x dim) and .g (k x 1) are an AND / polytope.
%   HalfSpace encodes  G * y <= g  (so each row is one linear "<= " constraint).
%
% Sign convention (from process_constraint):
%   (<= Y_a c)      -> row G(a)=+1,            g = c
%   (>= Y_a c)      -> row G(a)=-1,            g = -c
%   (<= Y_a Y_b)    -> row G(a)=+1, G(b)=-1,   g = 0
%   (>= Y_a Y_b)    -> row G(a)=-1, G(b)=+1,   g = 0
% (variable index a in Y_a is 0-based in the file -> column a+1 in G).
%
% To run: results = runtests('test_load_vnnlib_forms')

%% (a) input box + single (<= Y_0 0) -> 1 disjunct, 1 row
rng(42);
tf = [tempname '.vnnlib'];
fid = fopen(tf, 'w');
fprintf(fid, '(declare-const X_0 Real)\n(declare-const X_1 Real)\n');
fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n(declare-const Y_2 Real)\n');
fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
fprintf(fid, '(assert (<= X_1 2.0))\n(assert (>= X_1 -2.0))\n');
fprintf(fid, '(assert (<= Y_0 0.0))\n');
fclose(fid);
property = load_vnnlib(tf);
delete(tf);

% input box recovered exactly
assert(numel(property.lb) == 2, '(a) two input dims');
assert(isequal(double(property.lb(:)'), [-1 -2]), '(a) lb = [-1 -2]');
assert(isequal(double(property.ub(:)'), [ 1  2]), '(a) ub = [1 2]');
% one disjunct, one row
assert(numel(property.prop) == 1, '(a) single prop cell');
Hg = property.prop{1}.Hg;
assert(numel(Hg) == 1, '(a) exactly one HalfSpace (no OR)');
assert(size(Hg(1).G, 1) == 1, '(a) one constraint row');
assert(size(Hg(1).G, 2) == 3, '(a) three output dims wide');
% (<= Y_0 0) -> G(1)=+1, rest 0, g=0
assert(Hg(1).G(1, 1) == 1, '(a) Y_0 coeff = +1');
assert(all(Hg(1).G(1, 2:3) == 0), '(a) Y_1, Y_2 coeffs = 0');
assert(Hg(1).g(1) == 0, '(a) rhs g = 0');

%% (b) two standalone conjunct asserts (<= Y_0 0) and (>= Y_1 0) -> 1 disjunct, 2 rows
rng(42);
tf = [tempname '.vnnlib'];
fid = fopen(tf, 'w');
fprintf(fid, '(declare-const X_0 Real)\n');
fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n');
fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
fprintf(fid, '(assert (<= Y_0 0.0))\n');
fprintf(fid, '(assert (>= Y_1 0.0))\n');
fclose(fid);
property = load_vnnlib(tf);
delete(tf);

assert(numel(property.prop) == 1, '(b) single prop cell');
Hg = property.prop{1}.Hg;
assert(numel(Hg) == 1, '(b) one HalfSpace: two standalone asserts AND into one polytope');
assert(size(Hg(1).G, 1) == 2, '(b) two constraint rows (AND)');
% row 1: (<= Y_0 0) -> G=[+1 0], g=0
assert(isequal(Hg(1).G(1, :), [1 0]), '(b) row1 = [1 0] for Y_0<=0');
assert(Hg(1).g(1) == 0, '(b) row1 g = 0');
% row 2: (>= Y_1 0) -> G=[0 -1], g=-0=0
assert(isequal(Hg(1).G(2, :), [0 -1]), '(b) row2 = [0 -1] for Y_1>=0');
assert(Hg(1).g(2) == 0, '(b) row2 g = 0');

%% (c) pure (or (and >=)(and >=)(and >=)) -> 3 disjuncts, 1 row each
rng(42);
tf = [tempname '.vnnlib'];
fid = fopen(tf, 'w');
fprintf(fid, '(declare-const X_0 Real)\n');
fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n(declare-const Y_2 Real)\n');
fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
fprintf(fid, '(assert (or (and (>= Y_0 1.0)) (and (>= Y_1 1.0)) (and (>= Y_2 1.0))))\n');
fclose(fid);
property = load_vnnlib(tf);
delete(tf);

assert(numel(property.prop) == 1, '(c) single prop cell holding the OR');
Hg = property.prop{1}.Hg;
assert(numel(Hg) == 3, '(c) three OR disjuncts -> three HalfSpaces');
assert(size(Hg(1).G, 1) == 1, '(c) disjunct 1 has one row');
assert(size(Hg(2).G, 1) == 1, '(c) disjunct 2 has one row');
assert(size(Hg(3).G, 1) == 1, '(c) disjunct 3 has one row');
% (>= Y_k 1) -> G(k)=-1, g=-1
assert(isequal(Hg(1).G, [-1 0 0]) && Hg(1).g == -1, '(c) disjunct 1 = Y_0>=1');
assert(isequal(Hg(2).G, [0 -1 0]) && Hg(2).g == -1, '(c) disjunct 2 = Y_1>=1');
assert(isequal(Hg(3).G, [0 0 -1]) && Hg(3).g == -1, '(c) disjunct 3 = Y_2>=1');

%% (d) OR followed by TWO standalone conjuncts -> every disjunct gains BOTH rows
% This is the load_vnnlib bug-fix path generalized to >1 trailing conjunct.
rng(42);
tf = [tempname '.vnnlib'];
fid = fopen(tf, 'w');
fprintf(fid, '(declare-const X_0 Real)\n');
fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n(declare-const Y_2 Real)\n');
fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
fprintf(fid, '(assert (or (and (>= Y_0 1.0)) (and (>= Y_1 1.0))))\n');
fprintf(fid, '(assert (<= Y_2 5.0))\n');     % conjunct 1
fprintf(fid, '(assert (>= Y_2 -5.0))\n');    % conjunct 2
fclose(fid);
property = load_vnnlib(tf);
delete(tf);

Hg = property.prop{1}.Hg;
assert(numel(Hg) == 2, '(d) still two OR disjuncts (conjuncts must NOT add disjuncts)');
% each disjunct started with 1 row, grew by 2 -> 3 rows
assert(size(Hg(1).G, 1) == 3, '(d) disjunct 1 grew by 2 trailing conjunct rows');
assert(size(Hg(2).G, 1) == 3, '(d) disjunct 2 grew by 2 trailing conjunct rows');
% the original disjunct constraints (Y_0>=1 / Y_1>=1) still present in row 1
assert(isequal(Hg(1).G(1, :), [-1 0 0]) && Hg(1).g(1) == -1, '(d) disjunct1 row1 = Y_0>=1');
assert(isequal(Hg(2).G(1, :), [0 -1 0]) && Hg(2).g(1) == -1, '(d) disjunct2 row1 = Y_1>=1');
% both conjuncts on Y_2 appear in BOTH disjuncts (rows 2 and 3)
% (<= Y_2 5) -> [0 0 1], g=5 ;  (>= Y_2 -5) -> [0 0 -1], g=5
assert(isequal(Hg(1).G(2, :), [0 0 1]) && Hg(1).g(2) == 5,  '(d) disjunct1 conjunctA');
assert(isequal(Hg(1).G(3, :), [0 0 -1]) && Hg(1).g(3) == 5, '(d) disjunct1 conjunctB');
assert(isequal(Hg(2).G(2, :), [0 0 1]) && Hg(2).g(2) == 5,  '(d) disjunct2 conjunctA');
assert(isequal(Hg(2).G(3, :), [0 0 -1]) && Hg(2).g(3) == 5, '(d) disjunct2 conjunctB');
% the Y_2 conjunct variable appears in EVERY disjunct
assert(any(Hg(1).G(:, 3) ~= 0), '(d) Y_2 conjunct present in disjunct 1');
assert(any(Hg(2).G(:, 3) ~= 0), '(d) Y_2 conjunct present in disjunct 2');

%% (e) argmax-style (or (and (>= Y_5 Y_0)) ... ) -> N disjuncts, two-variable +1/-1 rows
% "Y_5 is NOT the strict argmin/argmax vs each other class" pattern.
rng(42);
tf = [tempname '.vnnlib'];
fid = fopen(tf, 'w');
fprintf(fid, '(declare-const X_0 Real)\n');
% Six output vars Y_0..Y_5 (unrolled: %%-section script tests cannot contain a for-loop).
fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n(declare-const Y_2 Real)\n');
fprintf(fid, '(declare-const Y_3 Real)\n(declare-const Y_4 Real)\n(declare-const Y_5 Real)\n');
fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
fprintf(fid, '(assert (or (and (>= Y_5 Y_0)) (and (>= Y_5 Y_1)) (and (>= Y_5 Y_2)) (and (>= Y_5 Y_3)) (and (>= Y_5 Y_4))))\n');
fclose(fid);
property = load_vnnlib(tf);
delete(tf);

Hg = property.prop{1}.Hg;
assert(numel(Hg) == 5, '(e) five OR disjuncts');
% (>= Y_5 Y_k) -> G(Y_5=col6) = -1, G(Y_k=col k+1) = +1, g = 0
% disjunct i compares Y_5 against Y_(i-1)
assert(size(Hg(1).G, 1) == 1, '(e) disjunct 1 single row');
assert(Hg(1).G(1, 6) == -1, '(e) disjunct1 Y_5 coeff = -1');
assert(Hg(1).G(1, 1) == 1,  '(e) disjunct1 Y_0 coeff = +1');
assert(Hg(1).g(1) == 0,     '(e) disjunct1 g = 0');
assert(Hg(2).G(1, 6) == -1 && Hg(2).G(1, 2) == 1, '(e) disjunct2 = Y_5>=Y_1');
assert(Hg(3).G(1, 6) == -1 && Hg(3).G(1, 3) == 1, '(e) disjunct3 = Y_5>=Y_2');
assert(Hg(4).G(1, 6) == -1 && Hg(4).G(1, 4) == 1, '(e) disjunct4 = Y_5>=Y_3');
assert(Hg(5).G(1, 6) == -1 && Hg(5).G(1, 5) == 1, '(e) disjunct5 = Y_5>=Y_4');
% each row is EXACTLY the two-variable +1/-1 pattern: one +1, one -1, rest 0
assert(sum(Hg(1).G(1, :) == 1) == 1 && sum(Hg(1).G(1, :) == -1) == 1, '(e) one +1 and one -1');
assert(nnz(Hg(1).G(1, :)) == 2, '(e) exactly two nonzero coeffs');

%% (f) single (<= Y_0 Y_1) -> 1 row; exact coefficients
rng(42);
tf = [tempname '.vnnlib'];
fid = fopen(tf, 'w');
fprintf(fid, '(declare-const X_0 Real)\n');
fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n');
fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
fprintf(fid, '(assert (<= Y_0 Y_1))\n');
fclose(fid);
property = load_vnnlib(tf);
delete(tf);

Hg = property.prop{1}.Hg;
assert(numel(Hg) == 1, '(f) single disjunct');
assert(size(Hg(1).G, 1) == 1, '(f) single row');
% (<= Y_0 Y_1) -> G(Y_0)=+1, G(Y_1)=-1, g=0
assert(Hg(1).G(1, 1) == 1,  '(f) Y_0 coeff = +1');
assert(Hg(1).G(1, 2) == -1, '(f) Y_1 coeff = -1');
assert(Hg(1).g(1) == 0,     '(f) g = 0');

%% Summary
disp('test_load_vnnlib_forms: all sections passed');
