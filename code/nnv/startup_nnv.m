fprintf('\nAdding dependencies to Matlab path...\n');

%% MATLAB Version Check
% NNV requires MATLAB R2023a (9.14) or newer
min_version = 'R2023a';
if verLessThan('matlab', '9.14')
    warning('NNV:MATLABVersion', ...
        'NNV requires MATLAB %s or newer. Current version: %s. Some features may not work.', ...
        min_version, version);
end

if ~exist("tbxmanager.m","file") % not added to the path
    addpath("tbxmanager");
end
% Restore toolbox paths. tbxmanager's enabled-toolbox list (tbxenabled.txt) is
% shared mutable state that is NOT concurrency-safe: parallel startup_nnv calls
% (e.g. batched/parallel verification jobs, as used when sweeping VNN-COMP
% benchmarks) or a process killed mid-write can corrupt it, after which
% 'restorepath' aborts with "Toolbox <name> is not installed". Self-heal by
% dropping any enabled entry whose toolbox is not actually installed, then retry,
% so a single corrupted line cannot take down every parallel worker's startup.
try
    tbxmanager restorepath
catch tbx_err
    warning('NNV:tbxRestore', ...
        'tbxmanager restorepath failed (%s); repairing tbxenabled.txt and retrying.', ...
        tbx_err.message);
    tbx_dir = fileparts(which('tbxmanager.m'));
    enabled_file = fullfile(tbx_dir, 'tbxenabled.txt');
    toolboxes_dir = fullfile(tbx_dir, 'toolboxes');
    if isfile(enabled_file)
        raw_lines = strsplit(strtrim(fileread(enabled_file)), newline);
        kept = {};
        for tbx_i = 1:numel(raw_lines)
            entry = strtrim(raw_lines{tbx_i});
            parts = strsplit(entry, ':');   % name:version:arch
            if numel(parts) == 3 && isfolder(fullfile(toolboxes_dir, parts{1}, parts{2}, parts{3}))
                kept{end+1} = entry; %#ok<AGROW>
            end
        end
        kept = unique(kept, 'stable');      % drop duplicates from a partial write
        fid = fopen(enabled_file, 'w');
        if fid > 0
            if ~isempty(kept)
                fprintf(fid, '%s\n', kept{:});
            end
            fclose(fid);
        end
    end
    tbxmanager restorepath                  % retry once on the repaired list
end

fprintf('\nAdding NNV to Matlab path...\n');
%mydir  = pwd;
%mydir
%idcs   = strfind(mydir,filesep);
%idcs
%newdir = mydir(1:idcs(end)-1);
%p = genpath(newdir); % generate a path that includes NNV folder and all folders below it
%addpath(p);
%cd(newdir);

addpath(pwd());
addpath(genpath(pwd()));
if is_codeocean()
    cd('/code/')
end

% display toolboxes to user
fprintf('\nNOTE: The following toolboxes are detected.\n');
ver()

%% Dependency Checking
% Check for required and optional toolboxes to help users debug issues
% rather than encountering cryptic errors when dependencies are missing

fprintf('\n--- NNV Dependency Check ---\n');

% Get installed toolbox names once
installed_toolboxes = {ver().Name};

% Critical dependencies - will likely cause significant errors if missing
% Format: {display_name, check_function}
% We use function existence as the primary check (most reliable)
critical_deps = {
    'Optimization Toolbox',                     'optimoptions'
    'Deep Learning Toolbox',                    'trainNetwork'
    'Image Processing Toolbox',                 'imresize'
    'Statistics and Machine Learning Toolbox',  'fitlm'
    'Parallel Computing Toolbox',               'parpool'
    'Computer Vision Toolbox',                  'detectSURFFeatures'
};

% Important but less critical - some features may not work
important_deps = {
    'Control System Toolbox',                   'ss'
    'Symbolic Math Toolbox',                    'sym'
    'System Identification Toolbox',            'iddata'
};

% Optional support packages - good to have for specific model imports
% These are support packages, not toolboxes, so we check function existence only
optional_deps = {
    'Deep Learning Toolbox Converter for ONNX Model Format',       'importONNXNetwork'
    'Deep Learning Toolbox Converter for PyTorch Models',          'importNetworkFromPyTorch'
    'Deep Learning Toolbox Converter for TensorFlow Models',       'importTensorFlowNetwork'
};

% Helper function to check toolbox/package availability
check_dependency = @(display_name, check_func) ...
    any(strcmp(installed_toolboxes, display_name)) || ...
    (exist(check_func, 'file') > 0);

% Check critical dependencies
missing_critical = {};
for i = 1:size(critical_deps, 1)
    if ~check_dependency(critical_deps{i,1}, critical_deps{i,2})
        missing_critical{end+1} = critical_deps{i,1};
    end
end

if ~isempty(missing_critical)
    fprintf('\n');
    warning('NNV:MissingCriticalDependency', ...
        ['CRITICAL: The following toolboxes are missing and NNV may not function correctly:\n' ...
         '  - %s\n' ...
         'Please install these toolboxes for full functionality.'], ...
        strjoin(missing_critical, '\n  - '));
end

% Check important dependencies
missing_important = {};
for i = 1:size(important_deps, 1)
    if ~check_dependency(important_deps{i,1}, important_deps{i,2})
        missing_important{end+1} = important_deps{i,1};
    end
end

if ~isempty(missing_important)
    fprintf('\n');
    warning('NNV:MissingImportantDependency', ...
        ['IMPORTANT: The following toolboxes are not installed:\n' ...
         '  - %s\n' ...
         'Some NNV features may be limited.'], ...
        strjoin(missing_important, '\n  - '));
end

% Check optional dependencies (info only, no warning)
missing_optional = {};
for i = 1:size(optional_deps, 1)
    if ~check_dependency(optional_deps{i,1}, optional_deps{i,2})
        missing_optional{end+1} = optional_deps{i,1};
    end
end

if ~isempty(missing_optional)
    fprintf('\nINFO: Optional support packages not installed (needed for specific model imports):\n');
    for i = 1:length(missing_optional)
        fprintf('  - %s\n', missing_optional{i});
    end
end

% Summary
n_critical = size(critical_deps, 1) - length(missing_critical);
n_important = size(important_deps, 1) - length(missing_important);
n_optional = size(optional_deps, 1) - length(missing_optional);

fprintf('\nDependency Summary: %d/%d critical, %d/%d important, %d/%d optional\n', ...
    n_critical, size(critical_deps, 1), ...
    n_important, size(important_deps, 1), ...
    n_optional, size(optional_deps, 1));
fprintf('--- End Dependency Check ---\n\n');

% import data structures from Hyst
javaaddpath(['engine', filesep, 'hyst', filesep, 'lib', filesep, 'Hyst.jar']);
import de.uni_freiburg.informatik.swt.spaceexboogieprinter.*;
import com.verivital.hyst.automaton.*;
import com.verivital.hyst.grammar.antlr.*;
import com.verivital.hyst.grammar.formula.*;
import com.verivital.hyst.importer.*;
import com.verivital.hyst.ir.*;
import com.verivital.hyst.junit.*;
import com.verivital.hyst.util.*;
import com.verivital.hyst.main.*;
import com.verivital.hyst.passes.*;
import com.verivital.hyst.printers.*;
import com.verivital.hyst.simulation.*;
import de.uni_freiburg.informatik.swt.sxhybridautomaton.*;
import de.uni_freiburg.informatik.swt.spaceexxmlprinter.*;
import de.uni_freiburg.informatik.swt.spaxeexxmlreader.*;

import com.verivital.hyst.automaton.*;
import com.verivital.hyst.grammar.antlr.*;
import com.verivital.hyst.grammar.formula.*;
import com.verivital.hyst.importer.*;
import com.verivital.hyst.ir.*;
import com.verivital.hyst.ir.base.*;
import com.verivital.hyst.ir.network.*;
import com.verivital.hyst.junit.*;
import com.verivital.hyst.main.*;
%import com.verivital.hyst.main.Hyst;
import com.verivital.hyst.outputparser.*;
import com.verivital.hyst.passes.*;
import com.verivital.hyst.passes.basic.*;
import com.verivital.hyst.passes.complex.*;
import com.verivital.hyst.passes.flatten.*;
%import com.verivital.hyst.passes.flatten.FlattenAutomatonPass;
import com.verivital.hyst.printers.*;
import com.verivital.hyst.python.*;
import com.verivital.hyst.simulation.*;
import com.verivital.hyst.util.*;

import de.uni_freiburg.informatik.swt.spaceexxmlprinter.*;
import de.uni_freiburg.informatik.swt.spaxeexxmlreader.*;
import de.uni_freiburg.informatik.swt.sxhybridautomaton.*;

%% Ready Confirmation
fprintf('\n');
fprintf('========================================\n');
fprintf('  %s Ready!\n', NNVVERSION());
fprintf('========================================\n');
fprintf('Get started: help Star, help NN\n');
fprintf('Tutorials:   examples/Tutorial/\n');
fprintf('Examples:    examples/\n');
fprintf('Support:     github.com/verivital/nnv/issues\n');
fprintf('\n');

%cd(mydir);
