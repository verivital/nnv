function ytrue = parse_argmax_vnnlib(vnnlibFile)
%PARSE_ARGMAX_VNNLIB Extract the predicted class index from an argmax-form VNNLIB.
%
%   Argmax-form output spec looks like:
%     (assert (or
%         (and (>= Y_0 Y_k))
%         (and (>= Y_1 Y_k))
%         ... (skipping (>= Y_k Y_k))
%         (and (>= Y_C Y_k))
%     ))
%   The right-hand `Y_k` (uniform across all disjuncts) is the predicted /
%   "true" class. Returns it as a 1-based MATLAB index.
%
%   Robustness: scans every (>= Y_<i> Y_<j>) match in the file, picks the
%   `j` that appears as the right-hand side most often (mode). For a
%   well-formed argmax spec this will be uniform and unambiguous.
%
%   Errors if no argmax pattern is detected (i.e. VNNLIB is half-space
%   form, not argmax form).

    if endsWith(vnnlibFile, ".gz")
        gunzip(vnnlibFile, tempdir);
        [~, base, ~] = fileparts(vnnlibFile);  % strips .gz
        vnnlibFile = fullfile(tempdir, base);
        cleanup = onCleanup(@() delete(vnnlibFile)); %#ok<NASGU>
    end

    fid = fopen(vnnlibFile, 'r');
    if fid == -1, error("parse_argmax_vnnlib: cannot open %s", vnnlibFile); end
    txt = fread(fid, Inf, 'char=>char')';
    fclose(fid);

    % Pattern: (>= Y_<i> Y_<j>) with optional whitespace.
    matches = regexp(txt, '\(\s*>=\s*Y_(\d+)\s+Y_(\d+)\s*\)', 'tokens');
    if isempty(matches)
        error("parse_argmax_vnnlib: no '(>= Y_i Y_j)' patterns found in %s — VNNLIB may be half-space form, not argmax", vnnlibFile);
    end

    rhs = cellfun(@(c) str2double(c{2}), matches);
    [u, ~, ic] = unique(rhs);
    counts = accumarray(ic, 1);
    [~, im] = max(counts);
    ytrue0 = u(im);          % 0-indexed class id (VNNLIB convention)
    ytrue  = ytrue0 + 1;     % 1-indexed MATLAB class id
end
