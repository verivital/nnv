function exportToC(obj, fname, dirname)
%MPT_EXPORTC Exports an explicit controller to C code
%
% ctrl.exportToC(obj, fname, dirname)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% The routine creates C-files for evaluation of the piecewise affine 
% control law. The first file constains a pure C script for evaluation of the
% explicit solution which can be ported to various platforms.
% The second file is an S-function implementation for verification of the
% controller under Simulink. 
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl   - MPT explicit controller
% fname  - Name of the file to be generated ("mpt_getInput" by default)
% dirname - Name of the directory that is generated ("mpt_explicit_controller"
%           by default)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
%     2012-2013 Revised by Martin Herceg, Automatic Control Laboratory, ETH
%          Zurich, herceg@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------



global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if nargin<2,
    fname = 'mpt_getInput';
else
    if isempty(fname)
        fname = 'mpt_getInput';
    end
    fname = strtrim(fname);
    if numel(fname)<=3
        error('The name of the file must have more than 3 characters.');
    end
    % get the short name if the full path is provided
    [~,fname] = fileparts(fname);   
    if ~isempty(regexp(fname,'\W', 'once'))
        error('The file name must contain only alphanumerical characters including underscore "_".');
    end    
end

if nargin<3
    dirname = 'mpt_explicit_controller';
else
    if isempty(dirname)
        dirname = 'mpt_explicit_controller';
    end
    dirname = strtrim(dirname);
    if numel(dirname)<1
        error('The name of the directory must have at least 1 character.');
    end
end

% check if there are multiple controllers
if numel(obj.feedback)>1
    error(['The code generation can be applied only for explicit controllers with a single ',...
        'feedback law. This controller has %d possible feedback laws. Use "PolyUnion/min" ',...
        'method to find single optimizer, i.e. \n\n',...
        '   ctrl.optimizer = ctrl.optimizer.min(''obj'')'],numel(obj.feedback));
end

% append the dirname to a current directory
dirname = [pwd, filesep, dirname];

% generate the files inside given directory
if ~mkdir(dirname)
    error('Could not create directory "%s".',dirname);
end
% create the file name with the path but without extension
fullfilename = [dirname, filesep, fname];

% call toC method of PolyUnion/BinTreePolyUnion
if isa(obj.optimizer,'BinTreePolyUnion')
    % no tie-breaking function for search tree
    obj.optimizer.toC('primal',fullfilename);
else
    % for all other cases use objective function as tie-breaking function
    obj.optimizer.toC('primal',fullfilename,'obj');
end
precision = MPTOPTIONS.modules.geometry.unions.(class(obj.optimizer)).toC.precision;
precision = strtrim(lower(precision));
if isempty(precision)
    precision = 'double';
elseif ~isequal(precision,'single') && ~isequal(precision,'double')
    error('The specified precision in the option can be either "single" or "double".');
end
if isequal(precision,'single')
    precision = 'float';
end


% sample time
Ts = obj.model.Ts;
if isempty(Ts)
    Ts = 1;
end

%% write mpt_getInput_sfunc.c
f2 = fopen([fullfilename,'_sfunc.c'], 'w');
if f2<0,
    error('Cannot open file for writing!');
end

part1={''
'/*'
'  Autogenerated C-code S-function for simulation of explicit controllers.'
'    '
'*/'
''
'/* Copyright (C) 2005 by Michal Kvasnica (michal.kvasnica@stuba.sk) '
'   Revised in 2012-2013 by Martin Herceg, Automatic Control Laboratory,'
'   ETH Zurich, herceg@control.ee.ethz.ch'
'*/'
''
'/*  This program is free software; you can redistribute it and/or modify'
'    it under the terms of the GNU General Public License as published by'
'    the Free Software Foundation; either version 2 of the License, or'
'    (at your option) any later version.'
''
'    This program is distributed in the hope that it will be useful,'
'    but WITHOUT ANY WARRANTY; without even the implied warranty of'
'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
'    GNU General Public License for more details.'
''
'    You should have received a copy of the GNU General Public License'
'    along with this program; if not, write to the Free Software'
'    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.'
'*/ '
''};

% write part1
for i=1:numel(part1)
    fprintf(f2,[part1{i},'\n']);
end

fprintf(f2,'/* Generated on %s by MPT %s */ \n\n',datestr(now), MPTOPTIONS.version);

% write the name of the file to include
fprintf(f2,'#define S_FUNCTION_NAME  %s_sfunc   /*Name of the S-function file*/\n',fname);
fprintf(f2,'#define S_FUNCTION_LEVEL 2		/*Level 2 S-functions allow multi-port status*/\n');
fprintf(f2,'#define MPT_TS %f		/* sample time */\n',Ts);
fprintf(f2,'#include "simstruc.h"			/*Header where different routines are defined */\n\n');
fprintf(f2, '#include "%s.c" /*inclusion of the evalution algorithm */\n\n',fname);


part2 = {''
''
'static void mdlInitializeSizes(SimStruct *S)'
'{'
'    ssSetNumSFcnParams(S, 0);'
''
'    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {'
'        return; /* Parameter mismatch will be reported by Simulink */'
'    }'
''
'    /* no states, all computation will be done in output section */'
'    ssSetNumContStates(S, 0);'
'    ssSetNumDiscStates(S, 0);'
'    '
'    /* one input port: state x(k) + vector of references */'
'    ssSetNumInputPorts(S, 1); '
'    '
'    /* one output port: control action */'
'    ssSetNumOutputPorts(S, 1);'
''
'    /* width of input vector - number of states of the original problem. '
'     * note that tracking includes additional states, here we do not consider them'
'     */'
'    ssSetInputPortWidth(S, 0, MPT_DOMAIN);'
'    '
'    ssSetInputPortDirectFeedThrough(S, 0, 1);'
'    '
'    /* width of output - number of control actions */'
'    ssSetOutputPortWidth(S, 0, MPT_RANGE);'
'    '
'    ssSetNumSampleTimes(S, 1);'
''
'    /* Take care when specifying exception free code - see sfuntmpl.doc */'
'    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);'
'    '
'}'
''
''
'static void mdlInitializeSampleTimes(SimStruct *S)'
'{'
'    /* set sampling time */'
'    ssSetSampleTime(S, 0, MPT_TS);'
'    ssSetOffsetTime(S, 0, 0.0);'
'}'
''
''
'#define MDL_INITIALIZE_CONDITIONS'
''
'static void mdlInitializeConditions(SimStruct *S)'
'{'
'}'
''
''
'static void mdlOutputs(SimStruct *S, int_T tid)'
'{'
'    real_T            	*u   = ssGetOutputPortRealSignal(S,0);'
'    InputRealPtrsType    Xin   = ssGetInputPortRealSignalPtrs(S,0);'
sprintf('    static %s U[MPT_RANGE], X[MPT_DOMAIN];',precision)
'    unsigned long int region;'
'    int i;'
''
'    /* prepare the input signal for passing to mpt_getInput function */'
'    for (i=0; i<MPT_DOMAIN; i++) {'
'         X[i] = *Xin[i];'
'    }'
''
''};

% write part2
for i=1:numel(part2)
    fprintf(f2,[part2{i},'\n']);
end

fprintf(f2,'    /* get control law */\n');
fprintf(f2,'    region = %s(X, U);\n\n',fname);

part3 = {
''
'    /* check if control law was found, if not, stop the simulation */'
'    if (region<1) {'
'        ssSetErrorStatus(S, "No feasible control law found!");'
'    }'
'    '
'    /* output control action */'
'    for (i=0; i<MPT_RANGE; i++) {'
'            u[i] = (real_T)U[i]; /* current output */'
'    }'
'}'
''
''
'/* Function: mdlTerminate ====================================================='
' * Abstract:'
' *    No termination needed, but we are required to have this routine.'
' */'
'static void mdlTerminate(SimStruct *S)'
'{'
'}'
''
'/*End of file necessary includes*/'
''
'#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */'
'#include "simulink.c"      /* MEX-file interface mechanism */'
'#else'
'#include "cg_sfun.h"       /* Code generation registration function */'
'#endif'
''};

% write part3
for i=1:numel(part3)
    fprintf(f2,[part3{i},'\n']);
end

% close the second file
fclose(f2);

end