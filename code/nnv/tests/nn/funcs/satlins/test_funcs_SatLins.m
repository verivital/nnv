%to run this as a test, use results_funcs_SatLins=runtests('test_funcs_SatLins')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%common variables
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star
I_poly = I.toPolyhedron;
sample_size=25;

%POTENTIAL QUESTION: do we want to draw different random variables for
%each test?

%POTENTIAL QUESTION: stepReach vs SatLinss, keep both? many of these
%tests seem identical except for the usage of a different method

%POTENTIAL QUSETION: need to figure out exactly what asserts to do for
%all of this


%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_Evaluate.m
%% test 1: SatLins Evaluate
x = [-1.5; 0.5; 2];
y = SatLins.evaluate(x);
assert(isequal(y, [-1; 0.5; 1]))


%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach.m
%% test 2:  SatLins Reach exact star
S = SatLins.reach(I, 'exact-star'); % exach reach set using star
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 3:  SatLins Reach approx star
S = SatLins.reach(I, 'approx-star'); % over-approximate reach set using star
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 4:  SatLins Reach abs dom
S = SatLins.reach(I, 'abs-dom'); % over-approximate reach set using abstract-domain
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);
%todo: figure our what to test here

%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach_abstract_domain.m
%% test 5: SatLins Abstract Domain
S = SatLins.reach_abstract_domain(I); % over-approximate reach set
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 6:  SatLins Reach
S = SatLins.reach(I); % exach reach set
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach_polyhedron_exact.m

%TEST DISABLED BECAUSE:
%it's a duplicate. see above (test 2)
% test 7: SatLins Reach Exact Star
%S = SatLins.reach(I, 'exact-star'); % exach reach set using star
%[res, fail_in, fail_out]=random_search(I, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%% test 7: SatLins Reach Exact Polyhedron
S = SatLins.reach(I_poly, 'exact-polyhedron');
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);
%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach_star_approx.m
%% test 8: SatLins Reach Star Approx
S = SatLins.reach_star_approx(I); % over-approximate reach set
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

% test 9: SatLins Reach
%TEST DISABLED BECAUSE:
%it's a duplicate. see above (test 6)
%S = SatLins.reach(I); % exach reach set
%[res, fail_in, fail_out]=random_search(I, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach_zono_approx.m
%% test 9: SatLins Reach Zono Approx
lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [2 1; 1.5 -2];
I_zono = I1.affineMap(A, []);

Z = SatLins.reach_zono_approx(I_zono); % over-approximation using zonotope is very conservative
I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
X = I2.sample(100);
Y = SatLins.evaluate(X);
S = SatLins.reach_star_approx(I2);
figure;
subplot(1, 2, 1);
I_zono.plot;
title('Input Set');
subplot(1, 2, 2);
Z.plot;
hold on;
S.plot;
hold on;
plot(Y(1, :), Y(2, :), '*b');
title('Output set, inside is Star, outside is Zonotope');


%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReach.m
%% test 10: SatLins Step Reach 
S = SatLins.stepReach(I, 1);
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReachAbstractDomain.m
%% test 11: SatLins Step Reach Abstract Domain
S = SatLins.stepReachAbstractDomain(I, 1);
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 12: SatLins Step Reach Abstract Domain
S = SatLins.stepReach(I,1);
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReachPolyhedronExact.m

%TEST DISABLED BECAUSE:
%it's a duplicate. see above (test 12)
% test 13: SatLins Step Reach
%S = SatLins.stepReach(I, 1);
%[res, fail_in, fail_out]=random_search(I, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReachStarApprox.m
%% test 13: SatLins Step Reach Star Approx
S = SatLins.stepReachStarApprox(I, 1);
[res, fail_in, fail_out]=random_search(I, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);


%TEST DISABLED BECAUSE:
%it's a duplicate. see above (test 12)
% test 14: SatLins Step Reach Star Approx
%S = SatLins.stepReach(I,1);
%[res, fail_in, fail_out]=random_search(I, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReachZonoApprox.m
%% test 14: SatLins Step Reach Zono Approx


lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [2 1; 1.5 -2];
I_zono = I1.affineMap(A, []);


figure;
I_zono.plot;
Z = SatLins.stepReachZonoApprox(I_zono, 1); % over-approximation using zonotope is very conservative
I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
S = SatLins.stepReachStarApprox(I2,1);
figure;
Z.plot;
hold on;
S.plot;



















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
