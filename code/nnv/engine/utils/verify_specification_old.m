function result = verify_specification_old(reachSet, property, reachOptions)
% verify a based on the interection between the reach set and halfspaces defining the property
%   (assumed to be the un_robust_region or unsafe region to prove)
% 
% Syntax:
%    [result] = verify_speficication(reachSet, property)
%
% Inputs:
%   - reachSet: computed output set of neural network (e.g. 1x1 Star)
%   - property: cell array of cell array defining all conditions that outputSet must satisfy
%
% Output: (should we change this to sat, unsat or unknown?)
%   - result: 0 ->  property failed
%             1 ->  property satisfied
%             2 ->  unknown
    
    arguments
        reachSet 
        property 
        reachOptions = []
    end

    R = reachSet;
    nr = length(R);    % number of output sets (for approx should be 1)
    
    % Process property to verify
    if iscell(property) % created from vnnlib (one or multiple halfSpaces)
        property = property{1};
        property = property.Hg; % property transformed into a HalfSpace(s)
    end

    % Begin verification
    np = length(property);
    if np == 1 % only one halfspace
        for k = 1:nr
            Set = R(k);
            if ~isa(Set, "Star")
                Set = Set.toStar;
            end
            S = Set.intersectHalfSpace(property.G, property.g); % compute intersection with unsafe/not robust region
            if isempty(S)
                result = 1; % no intersection with unsafe region = safe (unsat)
            else 
                result = 2; % intersection with safe and unsafe region = unknown or unsafe
                break;
            end
        end
    else
        for k = 1:nr
            if ~isa(R(k), "Star")
                R(k) = R(k).toStar;
            end
            if isa(R(k).V, 'gpuArray')
                R(k) = R(k).changeDevice('cpu');
            end
        end
        
        disp_stuff = ~isempty(reachOptions) && isfield(reachOptions, 'dis_opt') && strcmp(reachOptions.dis_opt, 'display');
        
        result = 1;
        fprintf("\n\n");
        for k = 1:nr % check every reach set vs OR property
            f = zeros(1, R(k).nVar, 'like', R(k).V); % objective function
            pred_lb = R(k).predicate_lb;
            pred_ub = R(k).predicate_ub;
            
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                max_calls = 1;
            else
                max_calls = poolobj.NumWorkers;
            end
            
            results = nan(np, 1);
            futures(1:np) = parallel.FevalFuture;
            num_unfetched_futures = 0;
            for cp = 1:np
                % if disp_stuff     % debug
                %     while length(futures(arrayfun(@(f) strcmp(f.State, "running"), futures))) > 0
                %         dum = 1;
                %     end
                %     fprintf("Processing cp: %d \t Memory used unnecessarily?: %.4g GB\n", cp, ((reachOptions.free_mem_B_before_verify_specification + reachOptions.free_swap_B_before_verify_specification) - (NN.get_free_mem_B + NN.get_free_swap_B))/2^30);
                % end
                if num_unfetched_futures >= max_calls
                    valid_futures = futures(arrayfun(@(f) ~strcmp(f.State, "unavailable"), futures));
                    [completedIdx, temp_result] = fetchNext(valid_futures);
                    results(completedIdx) = temp_result;
                    if disp_stuff   % debug
                        fprintf("\n\n%s \t ", datetime('now'));
                        fprintf("Index: %d \t Result: %d\n", completedIdx, temp_result);
                        num_running_futures = length(futures(arrayfun(@(f) strcmp(f.State, "running"), futures)));
                        fprintf("Running Futures: %d, ", num_running_futures);
                        free_mem_B = NN.get_free_mem_B;
                        fprintf("Free memory: %.4g GB", free_mem_B/2^30);
                        free_swap_B = NN.get_free_swap_B;
                        fprintf("Free swap: %.4g GB", free_swap_B/2^30);
                        mem_per_worker_right_now_B = (reachOptions.free_mem_B_before_verify_specification + reachOptions.free_swap_B_before_verify_specification - free_mem_B - free_swap_B)/num_running_futures;
                        fprintf("Memory per worker right now: %.4g GB\n", mem_per_worker_right_now_B/2^30);
                    end
                    if ~temp_result
                        result = 0;
                        cancel(valid_futures);
                        break;
                    end
                end
                
                
                if isfield(reachOptions, "free_mem_frac_for_verify_specification")
                    maxNumParWorkers = parcluster('local').NumWorkers;
                    if max_calls ~= maxNumParWorkers
                        pause_second = 0;
                        while NN.get_free_mem_frac < reachOptions.free_mem_frac_for_verify_specification || NN.get_idle_cpu < 0.1
                            if mod(pause_second, 60) == 0
                                fprintf("\nPausing due to free memory fraction being %.2g which is less than the specified threshold %.2g", NN.get_free_mem_frac, reachOptions.free_mem_frac_for_verify_specification);
                            end
                            pause(1);
                            pause_second = pause_second + 1;
                        end
                        if pause_second > 0
                            disp(' ')
                        end
                    end
                end
                
                [C_addition, d_addition] = Star.addition_to_C_d_by_intersection_with_halfspace(R(k), property(cp).G, property(cp).g);
                new_C = vertcat(R(k).C, C_addition);
                new_d = vertcat(R(k).d, d_addition);
                futures(cp) = parfeval(@Star.isEmptySet_Static, 1, f, new_C, new_d, pred_lb, pred_ub);
                num_unfetched_futures = num_unfetched_futures + 1;
            end
            
            if result
                num_unfinished_parfeval_calls = sum(isnan(results));
                for m = 1:num_unfinished_parfeval_calls
                    [completedIdx, temp_result] = fetchNext(futures);
                    results(completedIdx) = temp_result;
                    fprintf("Index: %d \t Result: %d\n", completedIdx, temp_result);
                    if ~temp_result
                        result = 0;
                        cancel(futures);
                        break;
                    end
                end
            end
            
            if ~result
                break;
            end
        end
        
        if ~result
            result = 2;
        end
    end

end % close function


