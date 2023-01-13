function [XLower, XUpper, label] = load_property(property)

% This will be used for both MATLAB and NNV. MATLAB can only verify
% robustess properties, so we can only use a subset of acas xu models
% Cannot use the verification function, will need to use the estimate
% bounds and then compare that to the property bounds

    switch property
        case 'prop1' % Safe if COC <= 1500 (assert (>= Y_0 3.991125645861615))
            XLower = [0.6; -0.5; -0.5; 0.45; -0.5];
            XUpper = [0.679857769; 0.5; 0.5; 0.5; -0.45];
            label = "prop1";
        case 'prop2' % Unsafe if COC is maximal
            XLower = [0.6; -0.5; -0.5; 0.45; -0.5];
            XUpper = [0.679857769; 0.5; 0.5; 0.5; -0.45];
            label = 1;
        case 'prop3' % Unsafe if COC is minimal
            XLower = [-0.303531156; -0.009549297; 0.493380324; 0.3; 0.3];
            XUpper = [-0.298552812; 0.009549297; 0.5; 0.5; 0.5];
            label = 'prop3';
        case 'prop4' % Unsafe if COC is minimal
            XLower = [-0.303531156; -0.009549297; 0.0; 0.318181818; 0.083333333];
            XUpper = [-0.298552812; 0.009549297; 0.0; 0.5; 0.166666667];
            label = 'prop4';
        otherwise
            error('Wrong property name. We are only evaluating prop2, prop3, and prop4');
    
    end

end

