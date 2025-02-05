function verOut = output_vs_mask(img, mask)

    diff_pixels = find(img ~= mask);

    verOut = img;

    for i=1:length(diff_pixels)
        if img(diff_pixels(i)) == 2
            continue
        elseif img(diff_pixels(i)) == 1
            verOut(diff_pixels(i)) = -1;
        else
            verOut(diff_pixels(i)) = -2;
        end
    end

end