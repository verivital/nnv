function python_path = cp_env()
%CP_ENV Returns the Python executable path for CP verification methods
%   Returns the path to the Python executable in the NNV virtual environment.
%   Supports Windows, macOS, and Linux platforms.
%
%   Usage:
%       python_path = cp_env()
%
%   Prerequisites:
%       1. Create virtual environment: python -m venv .venv
%       2. Install dependencies: pip install -r requirement.txt
%
%   See also: Prob_reach, verify_robustness_cp

    nnv_path = nnvroot();

    if ispc
        % Windows: .venv\Scripts\python.exe
        python_path = fullfile(nnv_path, '.venv', 'Scripts', 'python.exe');
    else
        % Unix (macOS/Linux): .venv/bin/python
        python_path = fullfile(nnv_path, '.venv', 'bin', 'python');
    end

    % Verify Python executable exists
    if ~isfile(python_path)
        error('cp_env:PythonNotFound', [...
            'Python virtual environment not found at: %s\n\n', ...
            'To set up CP verification, run these commands from the NNV root directory:\n', ...
            '  1. python -m venv .venv\n', ...
            '  2. %s\n', ...
            '  3. pip install -r requirement.txt\n\n', ...
            'See PYTHON_SETUP.md for detailed instructions.'], ...
            python_path, ...
            getActivateCommand());
    end
end

function cmd = getActivateCommand()
    if ispc
        cmd = '.venv\Scripts\activate';
    else
        cmd = 'source .venv/bin/activate';
    end
end
