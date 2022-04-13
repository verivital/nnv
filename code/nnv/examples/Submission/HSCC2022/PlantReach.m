function [stateSet,new_combos] = PlantReach(plant, init_set, combos, listofcommands)

% Compute state sets for the plant (all possible combs)
stateSet = [];
new_combos = [];
seen_combos = {};
combosFlag = false;
mmm = 1;

for j=1:size(combos,1)
        
    sFlag = false;
    
    % Look for repeated commands+branch to avoid doing the same operation
    % multiple times
    if combosFlag
        for k=1:size(seen_combos,1)
            if isequal(combos{j,1},seen_combos{k,1})
                if combos{j,4}==seen_combos{k,2}
                    sFlag = true;
                    break
                end
            end
        end
    end
        
    % Pass the command to the plant
    if ~sFlag
        pointSet = plant.stepReachStar(init_set(combos{j,4}),label2command(combos{j,1},listofcommands)); % PostProcessing(combos(j,3),listofcommands) = Up(j) it would simplify the code a little bit
        for jk=1:length(pointSet)
            for i=j:size(combos,1)
                if isequal(combos{i,1},combos{j,1})
                    if isequal(combos{i,4},combos{j,4})
                        ccc = {combos{i,3}, mmm, combos{i,3}, mmm};
                        new_combos = [new_combos; ccc];
                    end
                end
            end
            mmm = mmm + 1;
            
        end
        seen_combos = [seen_combos; {combos{j,1}, combos{j,4}}];
        combosFlag = true;
        stateSet = [stateSet pointSet];
    end
end
