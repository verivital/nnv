function score = robustScore(star_set,targetL)
    if length(star_set) > 1
        Xi = Star.merge_stars(star_set,1,'single');
    else
        Xi = star_set;
    end
    [mm,MM] = Xi.getRanges;
    [~, idx1] = max(mm);
    idxs = [];
    for idxm=1:length(mm)
        if MM(idxm) >= mm(idx1)
            idxs = [idxs idxm];
        end
    end
    if length(idxs) > 1
        if ismember(targetL,idxs)
            score = 2; % Unknown
        else
            score = 0; % Not robust
        end
    elseif ~isempty(idxs)
        if idxs == targetL
            score = 1; % Robust
        else
            score = 0; % Not Robust
        end
    else
        warning('No max index found... Probably an error estimating the ranges');
        score = 3;
    end
end

