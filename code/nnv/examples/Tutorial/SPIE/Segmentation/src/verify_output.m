function ver_img = verify_output(imageSet)

    % Get data ranges (output pixels)
    [lb,ub] = imageSet.getRanges;
    
    % 1) get correctly classified as 0 (background)
    ver_background = (ub <= 0);
    ver_background = ~ver_background;
    
    % 2) get correctly classified as 1 (lession)
    ver_lesion = (lb > 0);

    % 3) get verified output
    ver_img = 2*ones(size(lb)); % pixel = 2 -> unknown
    
    background = find(ver_background == 0); % 0
    ver_img(background) = 0;
    
    lesion = find(ver_lesion == 1); % 1
    ver_img(lesion) = 1;


end

