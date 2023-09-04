I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
S = Star(V', I.A, I.b); % star set 1

V = sample(S, 100);

figure;
Star.plot(S);
hold on;
plot(V(1, :), V(2,:), '*');


% sampling a star set 
function x = sample(obj, N)
    % @N: number of points in the samples
    % @x: a set of N sampled points in the star set 
    
    if N < 1
        error('Invalid number of samples');
    end
    
    % get predicate bounds
    lb = obj.predicate_lb;
    ub = obj.predicate_ub;
    if isempty(lb) && isempty(ub)
        [lb,ub] = obj.getRanges;
    end
    % input dimensions
    n = obj.dim;
    % random values (samples) within the predicate bounds
    values = (ub - lb).*rand(n,N) + lb;
    % generate random inputs from the initial star set
    x = obj.V(:,1) + obj.V(:,2:end)'*values;
    xT = obj.C*x <= obj.d;
    xT = all(xT == 1,1);
    idxs = find(xT == 1);
    x = x(:,idxs);
end