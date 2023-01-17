function spec = load_property(property,tool)

% This will be used for both MATLAB and NNV. MATLAB can only verify
% robustess properties, so we can only use a subset of acas xu models
% Cannot use the verification function, will need to use the estimate
% bounds and then compare that to the property bounds
switch tool 
    case "NNV"
        switch property
            case 'prop1' % Safe if COC <= 1500 (assert (>= Y_0 3.991125645861615))
                spec.XLower = [0.6; -0.5; -0.5; 0.45; -0.5];
                spec.XUpper = [0.679857769; 0.5; 0.5; 0.5; -0.45];
                spec.H = [1 0 0 0 0];
                spec.g = 3.991125645861615;
            case 'prop2' % Unsafe if COC is maximal
                spec.XLower = [0.6; -0.5; -0.5; 0.45; -0.5];
                spec.XUpper = [0.679857769; 0.5; 0.5; 0.5; -0.45];
                spec.H = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
                spec.g = [0;0;0;0];
            case 'prop3' % Unsafe if COC is minimal
                spec.XLower = [-0.303531156; -0.009549297; 0.493380324; 0.3; 0.3];
                spec.XUpper = [-0.298552812; 0.009549297; 0.5; 0.5; 0.5];
                spec.H = [-1 1 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1];
                spec.g = [0;0;0;0];
            case 'prop4' % Unsafe if COC is minimal
                spec.XLower = [-0.303531156; -0.009549297; 0.0; 0.318181818; 0.083333333];
                spec.XUpper = [-0.298552812; 0.009549297; 0.0; 0.5; 0.166666667];
                spec.H = [-1 1 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1];
                spec.g = [0;0;0;0];
            otherwise
                error('Wrong property name. We are only evaluating prop1, prop2, prop3, and prop4');
        end
    case "MATLAB"
        switch property
            case 'prop1' % Safe if COC <= 1500 (assert (>= Y_0 3.991125645861615))
                spec.XLower = [0.6; -0.5; -0.5; 0.45; -0.5];
                spec.XUpper = [0.679857769; 0.5; 0.5; 0.5; -0.45];
%                 spec.label = "prop1";
            case 'prop2' % Unsafe if COC is maximal
                spec.XLower = [0.6; -0.5; -0.5; 0.45; -0.5];
                spec.XUpper = [0.679857769; 0.5; 0.5; 0.5; -0.45];
                spec.label = 1;
            case 'prop3' % Unsafe if COC is minimal
                spec.XLower = [-0.303531156; -0.009549297; 0.493380324; 0.3; 0.3];
                spec.XUpper = [-0.298552812; 0.009549297; 0.5; 0.5; 0.5];
%                 spec.label = 'prop3';
            case 'prop4' % Unsafe if COC is minimal
                spec.XLower = [-0.303531156; -0.009549297; 0.0; 0.318181818; 0.083333333];
                spec.XUpper = [-0.298552812; 0.009549297; 0.0; 0.5; 0.166666667];
%                 spec.label = 'prop4';
            otherwise
                error('Wrong property name. We are only evaluating prop1, prop2, prop3, and prop4');
        end
    otherwise
        error("Only NNV and MATLAB allowed")
    
end

end

