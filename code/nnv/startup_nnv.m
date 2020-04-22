fprintf('\nAdding dependencies to Matlab path...\n');

tbxmanager restorepath

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
    
%cd(mydir);
