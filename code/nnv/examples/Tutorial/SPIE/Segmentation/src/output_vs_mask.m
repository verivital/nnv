function verOut = output_vs_mask(img, mask, pred)

    diff_pixels = find(img ~= mask);

    verOut = img;

    verOut(diff_pixels) = 2; % assume unknown to start

    for i=1:length(diff_pixels)
        if img(diff_pixels(i)) == 2
            continue
        elseif img(diff_pixels(i)) == 1 % false positive
            if pred(diff_pixels(i)) == 3
                verOut(diff_pixels(i)) = 3;
            end
            % verOut(diff_pixels(i)) = 2;
        else
            if pred(diff_pixels(i)) == 4
                verOut(diff_pixels(i)) = 4; % false negative
            end
            % verOut(diff_pixels(i)) = 2;
        end
    end

end