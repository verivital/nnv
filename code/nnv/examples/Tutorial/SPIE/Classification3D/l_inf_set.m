% Return an ImageStar of an linf attack
function I = l_inf_set(img, epsilon, max_value, min_value)
    imgSize = size(img);
    disturbance = epsilon * ones(imgSize, "like", img); % disturbance value
    lb = max(img - disturbance, min_value);
    ub = min(img + disturbance, max_value);
    I = VolumeStar(single(lb), single(ub)); % default: single (assume onnx input models)
end