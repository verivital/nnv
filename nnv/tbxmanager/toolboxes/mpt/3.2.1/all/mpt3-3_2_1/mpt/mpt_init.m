function mpt_init
% MPT3 installation script
%
% To install MPT3, simply call "mpt_init".
%
% The installation script automatically adds all required directories
% to your MATLAB path and looks up solvers installed on your machine.
%
% Note: after installation either save the updated path manually, or call
% mpt_init every time when MATLAB starts.

% global iteration counter (to detect how often mpt_solve is called)
global MPTOPTIONS

% remove MPT2 from the path
%
% we detect MPT2 path based on whether mpt_computeTrajectory can be seen,
% exploiting the fact that the file is in the base dir of MPT2 installation
mpt2_file_to_check = 'mpt_computeTrajectory';
if exist(mpt2_file_to_check, 'file')
	fprintf('Removing MPT2 from the path (we also need to clear your workspace)...\n');
	rmpath(genpath(fileparts(which(mpt2_file_to_check))));

	% we need to do "clear classes" to refresh the @polytope class.
	clear classes
	rehash toolbox
	% restore global variables
	global MPTOPTIONS
end

% check matlab compatibility. We need at least R2011a (7.12.0) or newer
% since we rely on the copyable class
try
	old_matlab = verLessThan('matlab', '7.12');
catch
	old_matlab = true;
end

% main path
mpt3_main_path = fileparts(which(mfilename));

%% check if YALMIP is installed
%p = which('yalmip');
%if isempty(p)
%    % download and install YALMIP
%    [filestr,status] = urlwrite('http://www.control.isy.liu.se/~johanl/YALMIP.zip','YALMIP.zip');
%    if status~=1
%        error('Could not download YALMIP from the internet, YALMIP cannot be installed automatically. Please, download and install YALMIP from "http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.Download" website');
%    end
%    % unzip the files
%    unzip(filestr,[mpt3_main_path,filesep,'modules',filesep,'yalmip']);
%    % delete the archive
%    if exist(filestr,'file')
%        delete(filestr);
%    end
%    % add to matlab path
%    ypath = [mpt3_main_path,filesep,'modules',filesep,'yalmip',filesep,'yalmip'];
%    fprintf('Adding YALMIP "%s" to the Matlab path.\n',ypath);
%    addpath(genpath(ypath));
%end

% add path to these subfolders if not set
subfolders = {'demos';
    'utils';
    'modules';
    'tests'};

% add current folder (in case we're somewhere else)
addpath(mpt3_main_path);

% get path settings
ps = path;
i=1;
while ~isempty(ps)
   [pstr{i}, ps] = strtok(ps, ':');
   i=i+1;
end
% add to path if not in there
for i=1:length(subfolders)
    str = [mpt3_main_path,filesep,subfolders{i}];
    if ~any(strcmp(str,pstr))
        addpath(genpath(str));
    end    
end

% create indexed access for matlab web browser
% bwpath=[fileparts(mpt3_main_path),filesep,'doc',filesep,'html'];
% if exist(bwpath,'dir')
%     builddocsearchdb(bwpath);
% end

% clear global variable MPTOPTIONS and functions which cache solvers
if ~isempty(MPTOPTIONS)
    delete(mptopt);
    clear global MPTOPTIONS
    clear mptopt
    clear mpt_subSolvers
end

% remove MPT2 settings. switch off warnings to prevent matlab from flooding
% us with messages that we are trying to redefine the @polytope class.
w = warning;
warning('off');
if ispref('MPT_toolbox')
	rmpref('MPT_toolbox');
end
warning(w);

% reset persistent preference settings
if ispref('MPT')
    rmpref('MPT');
end

% default settings are found in mptopt
Options = mptopt;


% prints the copyright notice
fprintf('\nMPT toolbox %s initialized...\n', Options.version);
fprintf('Copyright (C) 2003-%s by M. Kvasnica, C.N. Jones, and M. Herceg\n',datestr(now,'yyyy'));
fprintf('For news, visit the MPT web page at http://control.ee.ethz.ch/~mpt/\n');

disp(['            LP solver: ' Options.lpsolver]);
disp(['            QP solver: ' Options.qpsolver]);
if ~isempty(Options.solvers_list.MILP)
    disp(['          MILP solver: ' Options.milpsolver]);
end
if ~isempty(Options.solvers_list.MIQP)
    disp(['          MIQP solver: ' Options.miqpsolver]);
end
%disp(['Vertex enumeration: ', s_exsolver]);
disp([' parametric LP solver: ' Options.plpsolver]);
disp([' parametric QP solver: ' Options.pqpsolver]);
disp(' ');
disp('These default options can be changed. See "help mptopt" for more details.');

if old_matlab
	fprintf('\nMPT requires at least Matlab R2011a or newer for full operation.\n');
	fprintf('Some features (e.g. control synthesis) will not available on your system.\n\n');
end
