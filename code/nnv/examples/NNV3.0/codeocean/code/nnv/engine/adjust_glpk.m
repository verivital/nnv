function adjust_glpk()
    % ADJUST_GLPK - Fix glpkmex for compatibility with NNV
    %
    % This function modifies the glpk.m file to comment out a problematic
    % 'clear glpkcc' statement that can cause issues during reachability
    % analysis.
    %
    % Platform-aware: Works on Windows (win64), Linux (glnxa64), and
    % macOS (maci64/maca64).
    %
    % Usage:
    %   adjust_glpk()
    %
    % Note: This is typically called automatically by install.m

    try
        % Determine platform
        arch = computer('arch');  % 'glnxa64', 'win64', 'maci64', or 'maca64'

        % Construct platform-specific path
        glpk_base = fullfile('..', 'tbxmanager', 'toolboxes', 'glpkmex', '1.0', arch);
        glpk_folder = sprintf('glpkmex_1_0_%s', arch);
        filename = fullfile(glpk_base, glpk_folder, 'glpk.m');

        % Check if file exists
        if ~isfile(filename)
            fprintf('Note: glpkmex not found for platform %s\n', arch);
            return;
        end

        % Read the file
        fid = fopen(filename);
        if fid == -1
            fprintf('Warning: Could not open %s\n', filename);
            return;
        end
        cac = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '');
        fclose(fid);

        % Write modified file
        fid = fopen(filename, 'w');
        if fid == -1
            fprintf('Warning: Could not write to %s\n', filename);
            return;
        end

        change_here = 372;  % Line with 'clear glpkcc;'
        for jj = 1 : min(change_here-1, length(cac{1}))
            fprintf(fid, '%s\n', cac{1}{jj});
        end
        fprintf(fid, '%s\n', '%clear glpkcc;');  % Comment out the clear statement
        for jj = change_here+1 : length(cac{1})
            fprintf(fid, '%s\n', cac{1}{jj});
        end
        fclose(fid);

        fprintf('Successfully adjusted glpk.m for platform %s\n', arch);
    catch ME
        fprintf('Warning: Could not adjust glpk.m: %s\n', ME.message);
    end
end
