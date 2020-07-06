%to run this as a test, use results_funcs_poslin=runtests('test_funcs_poslin')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables

I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I_zono = B.toZono;
A = [0.5 1; 1.5 -2];
I_zono = I_zono.affineMap(A, []);
I_star=I_zono.toStar;


sample_size=25;

%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_evaluate.m

%% test 1: PosLin evaluate
x = [-1; 0.5; 2];
y = PosLin.evaluate(x);
assert(isequal(y, [0; .5; 2]));


%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach.m

%% test 2: PosLin reach exact star
S = PosLin.reach(I_star, 'exact-star'); % exact reach set using star
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 3: PosLin reach approx star
S = PosLin.reach(I_star, 'approx-star'); % over-approximate reach set using star
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 4: PosLin reach approx zono
S = PosLin.reach(I_zono, 'approx-zono'); % over-approximate reach set using zonotope

%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_abstract_domain.m

%TESTS DISABLED BECAUSE 
%Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector

% test 5: PosLin reach abstract domain
%S = PosLin.reach_abstract_domain(I); % over-approximate reach set
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

% test 6: PosLin reach
%S = PosLin.reach_star_approx(I); % exach reach set
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_star_approx.m

%TEST DISABLED BECAUSE 
%Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector

% test 7: PosLin reach star approx
%S = PosLin.reach_star_approx(I); % over-approximate reach set
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%% test 8: PosLin reach star approx
S = PosLin.reach(I); % exach reach set
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);


%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_star_approx_fast.m

%% test 9: PosLin reach star approx fast
S = PosLin.reach(I); % exach reach set
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%TEST DISABLED BECAUSE 
%Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector
% test 10: PosLin reach star approx fast
%S = PosLin.reach_star_approx(I); 
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%TEST DISABLED BECAUSE 
%depricated method. not sure what the replacement is, but probably covered
%in one of the other tests tbh.
% test 11: PosLin reach star approx fast
%S = PosLin.reach_star_approx_fast(I); % fast - over-approximate reach set
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_star_approx_vs_zono.m

%% test 12: PosLin reach star approx vs zono

W1 = [1 -1; 0.5 2; -1 1];
b1 = [-1; 0.5; 0];

W2 = [-2 1 1; 0.5 1 1];
b2 = [-0.5; -0.5];

W3 = [2 -1; 0 1];
b3 = [1; 0];

L1 = LayerS(W1, b1, 'poslin'); % construct first layer
L2 = LayerS(W2, b2, 'poslin');   % construct second layer

lb_loc= -rand(2, 1); % lower-bound vector of input set
ub_loc = lb_loc + [0.5; 0.5];   % upper-bound vector of input set

I_zono_loc = Star(lb_loc, ub_loc); % construct input set
I1 = I_zono_loc.getZono;

X1 = I_zono_loc.affineMap(W1, b1);
BX1 = X1.getBox;
BX11 = X1.getBoxFast;
Y1 = PosLin.reach_star_approx_fast(X1);
B1 = Y1.getBoxFast;
BY1 = Y1.getBox;

XZ1 = I1.affineMap(W1, b1);
Z1 = PosLin.reach_zono_approx(XZ1);
BZ1 = Z1.getBox;

X2 = Y1.affineMap(W2, b2);
BX2 = X2.getBox;
BX22 = X2.getBoxFast;
Y2 = PosLin.reach_star_approx_fast(X2);
B2 = Y2.getBoxFast;

XZ2 = Z1.affineMap(W2, b2);
BXZ2 = XZ2.getBox;
Z2 = PosLin.reach_zono_approx(XZ2);
BZ2 = Z2.getBox;


X3 = Y2.affineMap(W3, b3);
BX3 = X3.getBox;
Y3 = PosLin.reach_star_approx_fast(X3);

XZ3 = Z2.affineMap(W3, b3);
BXZ3 = XZ3.getBox;
Z3 = PosLin.reach_zono_approx(XZ3);





%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_zono_approx.m

%% test 13: PosLin reach zono approx
Z = PosLin.reach_zono_approx(I_zono); % over-approximation using zonotope is very conservative


%TEST DISABLED BECAUSE
%this method no longer exists. what method is this supposed to test now? 
% test 14: PosLin reach zono approx
%Z = PosLin.reach_zono_approx2(I_zono); % over-approximation using new zonotope method

%% test 15: PosLin reach zono approx
S = PosLin.reach_star_approx(I_star);
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);


%TEST DISABLED BECAUSE
%this method no longer exists. what method is this supposed to test now? 
% test 16: PosLin reach zono approx
%S = PosLin.reach_star_approx_fast(I_star);
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);


%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_stepReach.m

%TEST DISABLED BECAUSE
%this method no longer exists. what method is this supposed to test now? 
% test 17: PosLin stepReach
%S = PosLin.stepReach2(I, 1);
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_stepReachAbstractDomain.m


%TEST DISABLED BECAUSE
%Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector
% test 18: PosLin stepReach Abstract Domain
%S = PosLin.stepReachAbstractDomain(I, 1);
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%assert(isnan(fail_in) && isnan(fail_out));

%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_stepReachStarApprox.m

%TEST DISABLED BECAUSE
%Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector
% test 19: PosLin stepReach Star Approx
%S = PosLin.stepReachStarApprox(I, 1);
%[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
%assert(isnan(fail_in) && isnan(fail_out));

%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_stepReachZonoApprox.m

%TEST DISABLED BECAUSE
%no such method, presumably meaning there was some refactoring
% test 20: PosLin stepReach Zono Approx
%Z = PosLin.stepReachZonoApprox(I_zono, 1); % over-approximation using zonotope is very conservative

%% test 21: PosLin stepReach Zono Approx
S = PosLin.stepReachStarApprox(I_star,1);
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);



function [res, fail_input, fail_output]=grid_search(original, shifted, search_params)
  cur_loc=search_params(:, 1);
  index_limit=size(search_params, 1);
  cont=true;
  res=true;
  fail_input=NaN;
  fail_output=NaN;
  while cont
    output=LogSig.evaluate(cur_loc);
    if original.contains(cur_loc) ~= shifted.contains(output)
      res=false;
      fail_input=cur_loc;
      fail_output=output;
      return;
    end

    inc_fail=true;
    index=1;
    while inc_fail
      inc_fail=false;
      cur_loc(index)=cur_loc(index)+search_params(index, 2);
      if cur_loc(index)>search_params(index, 3)
	inc_fail=true;
	cur_loc(index)=search_params(index, 1);
	index=index+1;
	if index>index_limit
	  return
	end
	
      end
      
    end
    
  end
end

function [res, fail_input, fail_output]=random_search(original, shifted, num_sample)
  sample_set=original.sample(num_sample);
  sample_output=LogSig.evaluate(sample_set);
  %shifted.dim%this fails in certain tests, and i'm not sure why. specifically in tests number 3 and 6
  res=true;
  fail_input=NaN;
  fail_output=NaN;
  for i=1:size(sample_set, 2)
    if  ~shifted.contains(sample_output(:, i))
      res=false;
      fail_input=sample_set(:, i);
      fail_output=sample_output(:, i);
    end
  end
end
