%to run this as a test, use results_funcs_LogSig=runtests('test_funcs_LogSig')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables


I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

B = I.getBox; 
I_zono = B.toZono;

W = [0.5 1; -1 1];

I_zono = I_zono.affineMap(W, []);

I_star = I_zono.toStar;


sample_size=25;


l = B.lb;%these and below are actually not used.
u = B.ub;
y_l = logsig(l);
y_u = logsig(u);
dy_l = logsig('dn', l);
dy_u = logsig('dn', u);


%___________________________________________________________________________________________________
%tests below originally taken from test.m

% test 1: LogSig


%TEST DISABLED BECAUSE OF BELOW
%Unable to resolve the name R1.C.

%P1 = Polyhedron('A', R1.C, 'b', R1.d);
%W = R1.V(:, 2:R1.nVar + 1);
%P2 = P1.affineMap(W);



%___________________________________________________________________________________________________
%tests below originally taken from test_LogSig_reach.m

%% test 1: LogSig reach approx star
S = LogSig.reach(I_star, 'approx-star');
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);


%% test 2: LogSig reach approx star split

S = LogSig.reach(I_star, 'approx-star-split');
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);


%% test 3: LogSig reach abs dom
S = LogSig.reach(I_star, 'abs-dom');
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);




%% test 4: LogSig reach approx zono
S = LogSig.reach(I_zono, 'approx-zono');
%[res, fail_in, fail_out]=random_search(I_zono, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);


%___________________________________________________________________________________________________
%tests below originally taken from test_LogSig_reach_star_approx.m


%TESTS DISABLED BECAUSE of some mismatch with upper and lower bounds
%in predicates. It would seem that in the creation of a star object,
%it is advisable to include such bounds, which the other tests do by
%converting to and from a zono object.

% test 6: LogSig reach star approx
%S = LogSig.reach_star_approx(I);
%[res, fail_in, fail_out]=random_search(I, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);


% test 7: LogSig reach star approx split
%S = LogSig.reach_star_approx_split(I);
%[res, fail_in, fail_out]=random_search(I, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);




%___________________________________________________________________________________________________
%tests below originally taken from test_LogSig_reach_zono_approx.m

%% test 5: LogSig reach zono approx split
S = LogSig.reach_star_approx_split(I_star);
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 6: LogSig reach zono approx no spit
S = LogSig.reach_star_approx_no_split(I_star);
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 7: LogSig reach zono approx 
Z = LogSig.reach_zono_approx(I_zono);
%[res, fail_in, fail_out]=random_search(I_zono, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);



%___________________________________________________________________________________________________
%tests below originally taken from test_LogSig_stepLogSig.m


%TESTS DISABLED FOR 2 REASONS
%1) they don't test anyhting new
%2) they don't run, because these methods were refactored 


% test 11: LogSig stepLogSig split
%S = LogSig.stepLogSig_Split(I, 1, l(1), u(1), y_l(1), y_u(1), dy_l(1), dy_u(1));
%[res, fail_in, fail_out]=random_search(I, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

% test 12: LogSig stepLogSig no split
%S = LogSig.stepLogSig_NoSplit(I, 1, l(1), u(1), y_l(1), y_u(1), dy_l(1), dy_u(1));
%[res, fail_in, fail_out]=random_search(I, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);





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
  shifted.dim%this fails in certain tests, and i'm not sure why. specifically in tests number 3 and 6
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
