function [result,rT] = reach_collins(net, propertyFile)
%% Reachability analysis of collins benchmarks

% Load specification to verify
[lb_x, ub_x, lb_y, ub_y] = load_collins_vnnlib(propertyFile);
% Create reachability parameters and options
X = ImageStar(lb_x',ub_x');
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
% Compute reachability
rT = tic;
Y = net.reach(X, reachOptions); % Seems to be working
rT = toc(rT);
% Evaluate property
[y_lb,y_ub] = Y.getRanges;
disp(' ');
disp('===============================')
disp('RESULTS')
% disp(' ')
result = 1;
for k=1:length(y_lb)
    if ~((y_lb(k) >= lb_y(k) && y_ub(k) <= ub_y(k)))
        disp(' ');
        disp('Verification failed');
        result = 0;
    end
end
disp('Property is satisfied');
disp("Reachability computation time = "+string(rT) + " seconds")

end