function verOut = pred_vs_mask(img, mask)

    img = single(img);
    diff_pixels = find(img ~= mask);

    verOut = img;

    for i=1:length(diff_pixels)
        if img(diff_pixels(i)) == 2
            continue
        elseif img(diff_pixels(i)) == 1
            verOut(diff_pixels(i)) = 3; % false positive
            % verOut(diff_pixels(i)) = 2;
        else
            verOut(diff_pixels(i)) = 4; % false negative
            % verOut(diff_pixels(i)) = 2;
        end
    end

end