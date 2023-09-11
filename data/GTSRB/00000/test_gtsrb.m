N = 100;

l_files = ls;

for i = 1 : N
    fname = l_files(i,:);
    fname = strrep(fname, ' ', '');
    [fp,fname_no_ext,ext] = fileparts(fname);
    
    if ext == '.ppm'
        d = imread(fname);
        if size(d) == 0
            fname
            continue;
        end
        im(i).data = d;
    end
end