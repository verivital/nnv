    fid = fopen('../tbxmanager/toolboxes/glpkmex/1.0/glnxa64/glpkmex_1_0_glnxa64/glpk.m');
    cac = textscan( fid, '%s', 'Delimiter','\n','whitespace', '');
    fclose(fid)
    fid = fopen('../tbxmanager/toolboxes/glpkmex/1.0/glnxa64/glpkmex_1_0_glnxa64/glpk.m', 'w');
    change_here = 372;
    for jj = 1 : change_here-1
        fprintf(fid, '%s\n', cac{1}{jj});
    end
    fprintf(fid, '%s\n', '%clear glpkcc;');
    for jj = change_here+1: length(cac{1})
        fprintf(fid, '%s\n', cac{1}{jj});
    end
    fclose(fid);