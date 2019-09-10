function clearAllExceptVars(varNames)
    clearvars -except varNames;
    %clear GLOBAL
    %clear FUNCTIONS;
end