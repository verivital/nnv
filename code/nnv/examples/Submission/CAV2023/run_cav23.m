function run_cav23()
% Reproduce all CAV23 experiments

% Turn off figure display
set(0,'DefaultFigureVisible','off')

% Ensure folders exist before trying to save results

% 1) Neural ODEs
if ~exist('NeuralODEs/Comparison/nnvresults', 'dir')
   mkdir('NeuralODEs/Comparison/nnvresults')
end

if ~exist('NeuralODEs/RandomEx/results', 'dir')
   mkdir('NeuralODEs/RandomEx/results')
end

end

