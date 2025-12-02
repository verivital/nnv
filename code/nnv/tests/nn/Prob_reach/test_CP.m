% test_CP - Test Conformal Prediction verification
%
% This test requires:
%   1. Python virtual environment set up with torch, numpy, scipy
%   2. See PYTHON_SETUP.md for setup instructions

% Check if Python environment is configured for CP
try
    python_path = cp_env();
    python_exists = isfile(python_path);
catch
    python_exists = false;
end

if ~python_exists
    warning('test_CP:PythonNotConfigured', [...
        'Skipping test_CP: Python environment not configured.\n', ...
        'To enable CP verification, see PYTHON_SETUP.md for instructions.\n', ...
        'Quick setup:\n', ...
        '  1. python -m venv .venv\n', ...
        '  2. pip install -r requirement.txt']);
    % Create a passing assertion to mark test as skipped rather than failed
    assert(true, 'Test skipped - Python not configured');
else
    run("../../../examples/NN/cifar10/CP_verify_robustness.m")
end
