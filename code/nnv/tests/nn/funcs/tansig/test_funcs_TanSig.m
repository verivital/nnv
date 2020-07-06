%to run this as a test, use results_funcs_TanSig=runtests('test_funcs_TanSig')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.




I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star
B = I.getBox; 
I_zono = B.toZono;
sample_size=25;

W = [0.5 1; -1 1];

I_zono = I_zono.affineMap(W, []);

I_star = I_zono.toStar;

%__________________________________________________
%below taken from test_TanSig_reach.m
%% test 1: TanSig reach approx star
S = TanSig.reach(I_star, 'approx-star');
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);
%% test 2: TanSig reach abs dom
S = TanSig.reach(I_star, 'abs-dom');
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 3: TanSig reach approx zono
S = TanSig.reach(I_zono, 'approx-zono');
%[res, fail_in, fail_out]=random_search(I_zono, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);

%__________________________________________________
%below taken from test_TanSig_reach_star_zono_aprox.m
%% test 4: TanSig reach star aprox
S = TanSig.reach_star_approx(I_star);
[res, fail_in, fail_out]=random_search(I_star, S, sample_size);
if ~res
    display(fail_in)
    display(fail_out)
end
assert(res);

%% test 5: TanSig reach star aprox zono
S = TanSig.reach_zono_approx(I_zono);
%[res, fail_in, fail_out]=random_search(I_zono, S, sample_size);
%if ~res
%    display(fail_in)
%    display(fail_out)
%end
%assert(res);


%% test 6: TanSig reach star aprox relax
I = ExamplePoly.randVrep;
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

S = TanSig.reach_star_approx(I);
relaxFactor = 0.5;
S2 = TanSig.reach_star_approx(I, 'approx-star', relaxFactor);
X = I.sample(10);
Y = TanSig.evaluate(X);



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