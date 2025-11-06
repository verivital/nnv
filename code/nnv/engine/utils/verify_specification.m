function result = verify_specification(reachSet, property, reachOptions)
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
        R2 = Star;
        for k = 1:nr
            if ~isa(R(k), "Star")
                R2(k) = R(k).toStar;
            else
                R2(k) = R(k);
            end
        end
        R = R2;
        
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
        if disp_stuff
            fprintf("\n\n");
        end
        for k = 1:nr % check every reach set vs OR property
            f = zeros(1, R(k).nVar, 'like', R(k).V); % objective function
            pred_lb = R(k).predicate_lb;
            pred_ub = R(k).predicate_ub;
            C_const = parallel.pool.Constant(R(k).C);
            d_const = parallel.pool.Constant(R(k).d);
            pred_lb_const = parallel.pool.Constant(pred_lb);
            pred_ub_const = parallel.pool.Constant(pred_ub);
            f_const = parallel.pool.Constant(f);
            
            % poolobj = gcp('nocreate');
            % if isempty(poolobj)
            %     max_calls = 1;
            % else
            %     max_calls = poolobj.NumWorkers;
            % end
            
            results = nan(np, 1);
            futures = [];
            for cp = 1:np
                while ~isempty(futures) && NN.get_free_mem_frac < free_mem_frac_setting_hack
                    % disp("Pausing for 0.1 s because free memory " + NN.get_free_mem_B/2^30 + " is less than the desirable " + NN.get_total_mem_B*free_mem_frac_setting_hack/2^30)
                    pause(0.1);
                end
                % fetch some output(s)
                finished_futures_indices = find(arrayfun(@(f) strcmp(f.State, "finished"), futures));
                if ~isempty(finished_futures_indices)
                    [temp_results, finished_cp_nos] = fetchOutputs(futures(finished_futures_indices));
                    futures(finished_futures_indices) = [];
                    results(finished_cp_nos) = temp_results;
                    
                    if isfield(reachOptions, 'dis_opt') && strcmp(reachOptions.dis_opt, 'display')
                        fprintf("%s \t ", datetime('now'));
                        disp("Index: " + finished_cp_nos + ";    Result: " + temp_results);
                    end
                    
                    if ~all(temp_results)
                        result = 0;
                        cancel(futures);
                        break;
                    end
                end
                
                % if disp_stuff
                %     fprintf("Free memory before launching cp no. %d: %.2g GB", cp, NN.get_free_mem_B/2^30)
                % end
                
                % if isfield(reachOptions, "free_mem_frac_for_LPs")
                %     maxNumParWorkers = parcluster('local').NumWorkers;
                %     if max_calls ~= maxNumParWorkers
                %         pause_second = 0;
                %         while NN.get_free_mem_frac < reachOptions.free_mem_frac_for_LPs || NN.get_idle_cpu < 0.1
                %             if mod(pause_second, 60) == 0
                %                 fprintf("\nPausing due to free memory fraction being %.2g which is less than the specified threshold %.2g, or due to CPU usage being too high", NN.get_free_mem_frac, reachOptions.free_mem_frac_for_LPs);
                %             end
                %             pause(1);
                %             pause_second = pause_second + 1;
                %         end
                %         if pause_second > 0
                %             disp(' ')
                %         end
                %     end
                % end
                
                [C_addition, d_addition] = Star.addition_to_C_d_by_intersection_with_halfspace(R(k), property(cp).G, property(cp).g);
                
                new_code = 1;
                
                if ~new_code
                    new_C = vertcat(R(k).C, C_addition);
                    new_d = vertcat(R(k).d, d_addition);
                    future = parfeval(@Star.isEmptySet_Static, 2, f, new_C, new_d, pred_lb, pred_ub, cp);
                else
                    future = parfeval(@isEmptySet_Static_local, 2, f_const, C_const, d_const, pred_lb_const, pred_ub_const, cp, C_addition, d_addition);
                end
                
                if isempty(futures)
                    futures = future;
                else
                    futures(end + 1) = future;
                end
            end
            
            if result
                for m = 1:length(futures)
                    [completedIdx, temp_result, finished_cp_no] = fetchNext(futures);
                    results(finished_cp_no) = temp_result;
                    if isfield(reachOptions, 'dis_opt') && strcmp(reachOptions.dis_opt, 'display')
                        fprintf("%s \t ", datetime('now'));
                        fprintf("printing after parfeval launches finished; Index: %d \t Result: %d\n", finished_cp_no, temp_result);
                    end
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


function [bool, cp_no] = isEmptySet_Static_local(f_const, C_const, d_const, predicate_lb_const, predicate_ub_const, cp_no, C_addition, d_addition)
    [bool, cp_no] = Star.isEmptySet_Static(f_const.Value, [C_const.Value; C_addition], [d_const.Value; d_addition], predicate_lb_const.Value, predicate_ub_const.Value, cp_no);
end

