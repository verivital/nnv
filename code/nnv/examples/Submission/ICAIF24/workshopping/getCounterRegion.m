function counterExamples = getCounterRegion(inputSet, unsafeRegion, reachSet)
    % counterExamples = getCounterRegion(inputSet, unsafeRegion, reachSet)
    % NOTE: This is only to be used with exact-star method
    % unsafeRegion = HalfSpace (unsafe/undesired region)
    % inputSet = ImageStar/Star
    % reachSet = Star
    %
    % check the "safety" of the reachSet
    % Then, generate counterexamples

    % Initialize variables
    counterExamples = [];
    
    % Get halfspace variables
    G = unsafeRegion.G;
    g = unsafeRegion.g;
    
    % Check for valid inputs
    if ~isa(inputSet, "Star")
        error("Must be a Star");
    end
    if ~isa(reachSet, "Star")
        error("Must be Star or ImageStar");
    end

    % Begin counterexample computation
    n = length(reachSet); % number of stars in the output set
    V = inputSet.V;
    for i=1:n
        % Check for safety, if unsafe, add to counter
        if ~isempty(reachSet(i).intersectHalfSpace(G, g))
            counterExamples = [counterExamples Star(V, reachSet(i).C, reachSet(i).d,...
                reachSet(i).predicate_lb, reachSet(i).predicate_ub)];
        end
    end

end

