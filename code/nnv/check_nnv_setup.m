function check_nnv_setup()
    % CHECK_NNV_SETUP Quick diagnostic to identify NNV setup issues
    %
    % Usage:
    %   check_nnv_setup()
    %
    % This function performs a series of checks to verify that NNV is
    % properly installed and configured. It reports:
    %   - NNV and MATLAB versions
    %   - Path configuration status
    %   - Core class availability
    %   - Basic functionality tests
    %   - Available LP solvers
    %
    % Example:
    %   check_nnv_setup()  % Run all diagnostics
    %
    % See also: startup_nnv, NNVVERSION

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('        NNV Setup Diagnostic\n');
    fprintf('============================================\n\n');

    all_passed = true;

    %% Version Information
    fprintf('--- Version Information ---\n');

    % NNV Version
    try
        v = NNVVERSION();
        fprintf('[OK] NNV Version: %s\n', v);
    catch
        fprintf('[!!] NNV Version: Not available (NNVVERSION not found)\n');
        all_passed = false;
    end

    % MATLAB Version
    fprintf('[OK] MATLAB Version: %s\n', version);

    % Check MATLAB version requirement
    if verLessThan('matlab', '9.14')
        fprintf('[!!] WARNING: NNV requires MATLAB R2023a or newer\n');
        all_passed = false;
    else
        fprintf('[OK] MATLAB version meets requirements (R2023a+)\n');
    end
    fprintf('\n');

    %% Path Configuration
    fprintf('--- Path Configuration ---\n');

    % Check if NNV is on path
    if exist('Star', 'class') == 8
        fprintf('[OK] NNV classes are on path\n');
    else
        fprintf('[!!] NNV classes NOT on path - run startup_nnv first\n');
        all_passed = false;
    end

    % Check key directories
    nnv_root = fileparts(which('startup_nnv'));
    if ~isempty(nnv_root)
        fprintf('[OK] NNV root: %s\n', nnv_root);
    else
        fprintf('[!!] Cannot locate NNV root directory\n');
        all_passed = false;
    end
    fprintf('\n');

    %% Core Classes
    fprintf('--- Core Classes ---\n');
    core_classes = {'Star', 'ImageStar', 'Zono', 'Box', 'NN', 'NNCS'};
    for i = 1:length(core_classes)
        if exist(core_classes{i}, 'class') == 8
            fprintf('[OK] %s class available\n', core_classes{i});
        else
            fprintf('[!!] %s class NOT found\n', core_classes{i});
            all_passed = false;
        end
    end
    fprintf('\n');

    %% Basic Functionality Tests
    fprintf('--- Functionality Tests ---\n');

    % Test 1: Star set creation
    try
        S = Star([0;0], [1;1]);
        fprintf('[OK] Star set creation works\n');
    catch ME
        fprintf('[!!] Star set creation failed: %s\n', ME.message);
        all_passed = false;
    end

    % Test 2: Box to Star conversion
    try
        B = Box([0;0], [1;1]);
        S = B.toStar();
        fprintf('[OK] Box to Star conversion works\n');
    catch ME
        fprintf('[!!] Box to Star conversion failed: %s\n', ME.message);
        all_passed = false;
    end

    % Test 3: Simple layer evaluation
    try
        W = eye(2);
        b = [0; 0];
        layer = FullyConnectedLayer(W, b);
        out = layer.evaluate([1; 2]);
        if isequal(out, [1; 2])
            fprintf('[OK] FullyConnectedLayer evaluation works\n');
        else
            fprintf('[!!] FullyConnectedLayer evaluation gave unexpected result\n');
            all_passed = false;
        end
    catch ME
        fprintf('[!!] FullyConnectedLayer evaluation failed: %s\n', ME.message);
        all_passed = false;
    end

    % Test 4: ReLU layer
    try
        relu = ReluLayer();
        out = relu.evaluate([-1; 0; 1]);
        if isequal(out, [0; 0; 1])
            fprintf('[OK] ReluLayer evaluation works\n');
        else
            fprintf('[!!] ReluLayer evaluation gave unexpected result\n');
            all_passed = false;
        end
    catch ME
        fprintf('[!!] ReluLayer evaluation failed: %s\n', ME.message);
        all_passed = false;
    end
    fprintf('\n');

    %% LP Solvers
    fprintf('--- Available LP Solvers ---\n');

    % Check for linprog (Optimization Toolbox)
    if license('test', 'Optimization_Toolbox') && exist('linprog', 'file')
        fprintf('[OK] linprog (Optimization Toolbox)\n');
    else
        fprintf('[--] linprog not available\n');
    end

    % Check for GLPK
    if exist('glpk', 'file')
        fprintf('[OK] GLPK\n');
    else
        fprintf('[--] GLPK not available\n');
    end

    % Check for Gurobi
    if exist('gurobi', 'file')
        fprintf('[OK] Gurobi\n');
    else
        fprintf('[--] Gurobi not available\n');
    end
    fprintf('\n');

    %% Java/Hyst
    fprintf('--- Java/Hyst ---\n');
    try
        java_version = javaMethod('getProperty', 'java.lang.System', 'java.version');
        fprintf('[OK] Java runtime: %s\n', char(java_version));
    catch
        fprintf('[--] Java runtime not detected (Hyst features unavailable)\n');
    end
    fprintf('\n');

    %% Summary
    fprintf('============================================\n');
    if all_passed
        fprintf('  All checks PASSED - NNV is ready to use!\n');
    else
        fprintf('  Some checks FAILED - see above for details\n');
        fprintf('  Try running: startup_nnv\n');
    end
    fprintf('============================================\n');
    fprintf('\n');
    fprintf('Get started:\n');
    fprintf('  help Star      - Star set documentation\n');
    fprintf('  help NN        - Neural network documentation\n');
    fprintf('  examples/      - Example scripts\n');
    fprintf('\n');
end
