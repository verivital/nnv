function display(P)
%
% display function
%

if numel(P)==0
    fprintf('Empty polyhedron array.\n');
    return;
elseif numel(P) > 1
    fprintf('Array of %i polyhedra.\n', numel(P));
    return
end

if isempty(P.H_int) && isempty(P.V_int) && isempty(P.R_int)
    if size(P.He_int,1) == 0
        fprintf('Empty polyhedron in R^%i\n', P.Dim);
    else
        fprintf('Affine set with %i equations in R^%i\n', size(P.He,1), P.Dim);
    end
    return
end
fprintf('Polyhedron in R^%i with representations:\n', P.Dim);
fprintf('    H-rep ');
if P.hasHRep
    if P.irredundantHRep, fprintf('%-13s', '(irredundant)');
    else                  fprintf('%-13s', '(redundant)');
    end
    fprintf(' : Inequalities %3i | Equalities %3i\n', size(P.H,1), size(P.He,1));
else
    fprintf('%-13s : Unknown (call computeHRep() to compute)\n', '');
end
fprintf('    V-rep ');
if P.hasVRep
    if P.irredundantVRep, fprintf('%-13s', '(irredundant)');
    else                  fprintf('%-13s', '(redundant)');
    end
    fprintf(' : Vertices %3i | Rays %3i\n', size(P.V,1), size(P.R,1));
else
    fprintf('%-13s : Unknown (call computeVRep() to compute)\n', '');
end

% display attached functions (implemented in ConvexSet/displayFunctions)
P.displayFunctions;

end
